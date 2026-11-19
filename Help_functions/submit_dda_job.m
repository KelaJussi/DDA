function entry = submit_dda_job(Input, Array, fullpath, filename, kind)
% SUBMIT_DDA_JOB  Start a DDA run without waiting for it, and register it.
%
%   Input.parallel = "cluster"     batch job via MATLAB Parallel Server (Input.cluster.*).
%                                  Setting Input.cluster.name = "Processes" runs the same
%                                  kind of job on this computer, with its own pool.
%   Input.parallel = "background"  parfeval on a local pool (Input.background.*); each
%                                  run uses one worker, several runs go in parallel.
%
% Returns the registry entry. Fetch finished jobs with collect_dda_jobs, see
% them with dda_jobs.

entry = struct( ...
    'key',        string(datetime('now','Format','yyyyMMdd_HHmmss')) + "_" + sprintf('%04d', randi(9999)), ...
    'type',       "", ...
    'profile',    "", ...
    'job_id',     NaN, ...
    'delete_job', true, ...
    'submitted',  datetime('now'), ...
    'fullpath',   string(fullpath), ...
    'filename',   string(filename), ...
    'kind',       string(kind), ...
    'plot_on_collect', getopt(getopt(Input,'jobs',struct()), 'plot_on_collect', false), ...
    'request',    struct('wall_time',"", 'num_workers',1, 'num_threads',NaN, 'partition',""), ...
    'Input',      Input);

dda_size_check(Input, Array);       % stop here, with a reason, if the run cannot fit

switch string(Input.parallel)
    case "cluster"
        cc = struct_defaults(getopt(Input,'cluster',struct()), cluster_defaults());
        cl = parcluster(char(cc.name));
        if ~isempty(cc.num_threads), try, cl.NumThreads = cc.num_threads; catch, end, end
        apply_additional_properties(cl, cc);
        if cc.save_profile, saveProfile(cl); end

        % A local profile ("Processes") shares this computer's files: give its
        % workers this MATLAB's path instead of attaching files.
        add_path = cc.auto_add_client_path || ~isprop(cl, 'AdditionalProperties');
        args = {'Pool', cc.num_workers, ...
                'CurrentFolder', char(cc.current_folder), ...
                'AutoAddClientPath', add_path, ...
                'AutoAttachFiles', cc.auto_attach_files, ...
                'CaptureDiary', true};
        if ~isempty(cc.attached_files),   args = [args, {'AttachedFiles',   cellstr(cc.attached_files)}];   end
        if ~isempty(cc.additional_paths), args = [args, {'AdditionalPaths', cellstr(cc.additional_paths)}]; end

        % "threads": inside the job, the eigen solver uses the job's own pool
        try
            job = batch(cl, @dda_solve_job, 1, {Input, Array, "threads", true}, args{:});
        catch err
            if ~is_connection_error(err), rethrow(err); end
            % A dead SSH connection (e.g. after the certificate expired): renew it
            % and submit once more, instead of restarting MATLAB.
            fprintf(['The connection to the cluster was lost (expired SSH certificate?). ' ...
                     'Reconnecting and submitting again...\n']);
            if ~dda_cluster_reconnect(cc.name)
                error('submit_dda_job:connection', ['Could not renew the connection to "%s". Renew the ' ...
                    'SSH certificate and run dda_cluster_reconnect, or restart MATLAB.\nOriginal error: %s'], ...
                    cc.name, err.message);
            end
            cl = parcluster(char(cc.name));
            if ~isempty(cc.num_threads), try, cl.NumThreads = cc.num_threads; catch, end, end
            apply_additional_properties(cl, cc);
            job = batch(cl, @dda_solve_job, 1, {Input, Array, "threads", true}, args{:});
        end
        entry.type = "cluster";  entry.profile = string(cc.name);
        entry.job_id = job.ID;   entry.delete_job = cc.delete_job;
        entry.request = struct('wall_time', string(cc.wall_time), 'num_workers', cc.num_workers, ...
                               'num_threads', cl.NumThreads, 'partition', string(cc.partition));
        if ~isprop(cl, 'AdditionalProperties'), entry.request.wall_time = ""; end   % local: no limit
        dda_job_registry("add", entry);
        fprintf('Submitted job %d to "%s" (%d workers + 1, %d threads each, wall time %s):\n  %s\n', ...
            job.ID, cc.name, cc.num_workers, cl.NumThreads, string(cc.wall_time), filename);

    case "background"
        bo = struct_defaults(getopt(Input,'background',struct()), ...
                             struct('pool', "Processes", 'num_workers', []));
        pool = gcp('nocreate');
        if isempty(pool)
            if isempty(bo.num_workers), pool = parpool(char(bo.pool));
            else,                       pool = parpool(char(bo.pool), bo.num_workers);
            end
        elseif isa(pool, 'parallel.ThreadPool') ~= (lower(string(bo.pool)) == "threads")
            % MATLAB allows one pool at a time; an open pool (e.g. kept by a
            % "threads" run) is used as it is.
            fprintf(['Note: using the open %s instead of Input.background.pool = "%s". ' ...
                     'To switch: delete(gcp) and submit again.\n'], class(pool), bo.pool);
        end
        F = parfeval(pool, @dda_solve_job, 1, Input, Array, "none", false);
        entry.type = "background";  entry.profile = string(class(pool));
        dda_job_registry("add", entry, F);
        fprintf('Started background job on the local pool (%d workers):\n  %s\n', pool.NumWorkers, filename);

    otherwise
        error('submit_dda_job:mode', 'Input.parallel must be "cluster" or "background" to submit a job.');
end
end

function tf = is_connection_error(err)
% A lost SSH connection to the cluster, anywhere in the chain of causes
msgs = strings(0);
todo = {err};
while ~isempty(todo)
    e = todo{1};  todo(1) = [];
    msgs(end+1) = string(e.message); %#ok<AGROW>
    todo = [todo, e.cause(:).']; %#ok<AGROW>
end
txt = lower(strjoin(msgs, " "));
tf = any(contains(txt, ["channel-open", "failed to start mirror", ...
         "failed to execute command on remote host", "connection is not established", ...
         "unable to connect", "connection reset", "broken pipe"]));
end

%% ========================================================================
%  Cluster options (moved here from Calculate_DDA_eigenmodes)
%% ========================================================================
function d = cluster_defaults()
% Every field can be overridden through Input.cluster in DDA_Master.
d = struct( ...
    'name',                 "Roihu R2026a", ... % parcluster profile ("Processes" = this computer)
    'wait',                 false, ...          % true: block until done, then save and plot
    'wall_time',            "0:30:00", ...
    'num_threads',          4, ...
    'num_workers',          10, ...
    'partition',            "", ...
    'project',              "", ...
    'mem_per_cpu',          "", ...
    'email',                "", ...
    'submit_args',          "", ...
    'additional_properties',struct(), ...
    'property_names',       struct(), ...
    'current_folder',       ".", ...
    'auto_add_client_path', false, ...
    'auto_attach_files',    true, ...
    'attached_files',       {{}}, ...
    'additional_paths',     {{}}, ...
    'delete_job',           true, ...
    'save_profile',         false);
end

function apply_additional_properties(cl, cc)
% Maps the frontend names to the plugin's AdditionalProperties names (the
% MathWorks Slurm plugin used on Roihu; override via Input.cluster.property_names,
% e.g. project -> 'ComputingProject' for the older CSC profiles).
% Local profiles ("Processes") have no AdditionalProperties: nothing to set.
if ~isprop(cl, 'AdditionalProperties'), return; end
names = struct('wall_time','WallTime', 'partition','Partition', ...
               'project','AccountName', 'mem_per_cpu','MemPerCPU', ...
               'email','EmailAddress', 'submit_args','AdditionalSubmitArgs');
names = struct_defaults(cc.property_names, names);
keys = fieldnames(names);
for n = 1:numel(keys)
    val = cc.(keys{n});
    if ~is_unset(val), set_cluster_property(cl, names.(keys{n}), val); end
end
extra = fieldnames(cc.additional_properties);
for n = 1:numel(extra)
    set_cluster_property(cl, extra{n}, cc.additional_properties.(extra{n}));
end
end

function set_cluster_property(cl, name, val)
known = {};
try, known = fieldnames(cl.AdditionalProperties); catch
    try, known = properties(cl.AdditionalProperties); catch, end
end
if ~isempty(known) && ~ismember(name, known)
    warning('submit_dda_job:unknownClusterProperty', ...
        ['"%s" is not an existing AdditionalProperties entry of this profile; setting it anyway. ' ...
         'Map the right name with Input.cluster.property_names.'], name);
end
if isstring(val), val = char(val); end
cl.AdditionalProperties.(name) = val;
end

function tf = is_unset(v)
tf = isempty(v) || ((isstring(v) || ischar(v)) && strlength(string(v)) == 0);
end

function s = struct_defaults(s, d)
f = fieldnames(d);
for n = 1:numel(f)
    if ~isfield(s,f{n}), s.(f{n}) = d.(f{n}); end
end
end

function v = getopt(s,name,default)
if isstruct(s) && isfield(s,name), v = s.(name); else, v = default; end
end
