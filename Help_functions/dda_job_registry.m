function varargout = dda_job_registry(action, varargin)
% DDA_JOB_REGISTRY  Book-keeping of DDA jobs that run without blocking MATLAB.
% Used by submit_dda_job, collect_dda_jobs, dda_jobs and DDA_calc_looped;
% you normally do not call it yourself.
%
% Each job is one small file <key>.mat in a local folder that OneDrive does
% not sync (one file per job, so two MATLAB sessions or a sync client can
% never produce "conflicted copies" of a shared index):
%   Windows: %LOCALAPPDATA%\DDA_jobs      other: ~/.local/share/DDA_jobs
% Cluster jobs can be found again after MATLAB restarts (findJob). Background
% (parfeval) jobs live in this MATLAB session only: their Future objects are
% kept by dda_job_futures (locked against "clear all"); after a restart, or
% when the pool is deleted, such jobs are reported as lost. This function is
% deliberately not locked, so edits to it take effect immediately.
%
% Actions
%   folder                     path of the registry folder
%   add,    entry, future      register a job (future only for background jobs)
%   list                       all entries, oldest first (struct array)
%   find,   fullpath           the entry of a result file, or []
%   state,  entry              "queued" "running" "finished" "failed" "timeout" "lost" "missing"
%                              "timeout": MATLAB still says running, but the job has run
%                              past its wall time (+2 min): Slurm stopped it without
%                              MATLAB noticing.
%   wait,   entry              block until the job has finished or failed
%   fetch,  entry              the job's output (struct R from dda_solve_job)
%   message, entry             error message and the end of the job's diary
%   times,  entry              struct with submitted, started, finished (datetime,
%                              NaT when unknown) as recorded by MATLAB
%   local_bytes, entry         bytes of the job's files in the local job storage
%                              (what was downloaded from the cluster); NaN if unknown
%   cancel, entry              stop the job (it stays registered)
%   remove, entry              forget the job; deletes a cluster job if entry.delete_job

persistent clusters
if isempty(clusters), clusters = containers.Map('KeyType','char','ValueType','any'); end

switch lower(string(action))
    case "folder"
        varargout{1} = registry_folder();

    case "add"
        entry = varargin{1};
        save(entry_file(entry.key), 'entry');
        if numel(varargin) > 1 && ~isempty(varargin{2}), dda_job_futures('set', entry.key, varargin{2}); end

    case "list"
        d = dir(fullfile(registry_folder(), '*.mat'));
        entries = [];
        for n = 1:numel(d)
            try
                S = load(fullfile(d(n).folder, d(n).name), 'entry');
                entries = [entries, normalize_entry(S.entry)]; %#ok<AGROW>
            catch
                warning('dda_job_registry:unreadable', 'Skipping unreadable registry file %s', d(n).name);
            end
        end
        if ~isempty(entries)
            [~, o] = sort([entries.submitted]); entries = entries(o);
        end
        varargout{1} = entries;

    case "find"
        target = string(varargin{1});
        entries = dda_job_registry("list");
        varargout{1} = [];
        for n = 1:numel(entries)
            if string(entries(n).fullpath) == target, varargout{1} = entries(n); return; end
        end

    case "state"
        varargout{1} = job_state(varargin{1});

    case "wait"
        h = job_handle(varargin{1});
        if ~isempty(h), wait(h); end

    case "local_bytes"
        e = varargin{1};  varargout{1} = NaN;
        if e.type == "cluster"
            try
                p = char(e.profile);
                if ~isKey(clusters, p), clusters(p) = parcluster(p); end
                loc = clusters(p).JobStorageLocation;
                if isstruct(loc)                     % {windows, unix} form
                    if ispc, loc = loc.windows; else, loc = loc.unix; end
                end
                varargout{1} = folder_bytes(char(loc), sprintf('Job%d', e.job_id));
            catch
            end
        end

    case "fetch"
        e = varargin{1};  h = job_handle(e);
        if isempty(h), error('dda_job_registry:gone', 'Job %s no longer exists.', e.key); end
        if e.type == "cluster"
            out = fetchOutputs(h);  varargout{1} = out{1};
        else
            varargout{1} = fetchOutputs(h);
        end

    case "message"
        e = varargin{1};  h = job_handle(e);  msg = "";  dia = "";
        if ~isempty(h)
            if e.type == "cluster"
                t = h.Tasks(1);
                if ~isempty(t.Error), msg = string(t.Error.message); end
                dia = string(t.Diary);
                if msg == "" && strcmp(h.State,'failed')
                    msg = "the scheduler reports the job as failed (wall time or memory exceeded?)";
                end
            else
                if ~isempty(h.Error), msg = string(h.Error.message); end
                dia = string(h.Diary);
            end
        end
        lines = splitlines(strtrim(dia));
        varargout{1} = msg;
        varargout{2} = strjoin(lines(max(1,end-9):end), newline);   % last 10 diary lines

    case "times"
        h = job_handle(varargin{1});
        t = struct('submitted', NaT, 'started', NaT, 'finished', NaT);
        if ~isempty(h)
            t.submitted = job_time(h, ["SubmitDateTime","CreateDateTime"]);
            t.started   = job_time(h, "StartDateTime");
            t.finished  = job_time(h, "FinishDateTime");
        end
        varargout{1} = t;

    case "cancel"
        h = job_handle(varargin{1});
        if ~isempty(h), cancel(h); end

    case "remove"
        e = varargin{1};
        if e.type == "cluster" && e.delete_job
            h = job_handle(e);
            if ~isempty(h), try, delete(h); catch, end, end
        end
        dda_job_futures('remove', e.key);
        f = entry_file(e.key);
        if isfile(f), delete(f); end

    otherwise
        error('dda_job_registry:action', 'Unknown action "%s".', action);
end

    %% --------------------------------------------------------------------
    function h = job_handle(e)
        h = [];
        if e.type == "cluster"
            p = char(e.profile);
            if ~isKey(clusters, p), clusters(p) = parcluster(p); end
            h = findJob(clusters(p), 'ID', e.job_id);
        else
            h = dda_job_futures('get', e.key);
        end
    end

    function s = job_state(e)
        h = job_handle(e);
        if isempty(h)
            if e.type == "cluster", s = "missing"; else, s = "lost"; end
            return
        end
        switch string(h.State)
            case {"pending","queued"}, s = "queued";
            case "running"
                s = "running";
                if e.type == "cluster" && past_wall_time(e, h), s = "timeout"; end
            case "finished"
                if e.type == "cluster", err = h.Tasks(1).Error; else, err = h.Error; end
                if isempty(err), s = "finished"; else, s = "failed"; end
            case "failed",             s = "failed";
            otherwise,                 s = "lost";        % 'unavailable', deleted pool
        end
    end
end

%% ------------------------------------------------------------------------
function tf = past_wall_time(e, h)
% A "running" cluster job that has run longer than its wall time plus a margin
% was stopped by Slurm: the killed worker cannot tell MATLAB, so the job state
% stays "running". The margin covers Slurm's grace period before killing and
% clock differences between the cluster and this computer.
tf = false;
req = NaN;
if isfield(e, 'request') && isfield(e.request, 'wall_time'), req = dda_walltime_seconds(e.request.wall_time); end
if isnan(req), return; end
started = job_time(h, "StartDateTime");
if isnat(started), return; end
elapsed = seconds(datetime('now', 'TimeZone', started.TimeZone) - started);
tf = elapsed > req + 120;
end

function e = normalize_entry(e)
% Entries saved by older versions lack newer fields; add them so that all
% entries can be put in one struct array.
if ~isfield(e, 'request')
    w = "";
    if isfield(e,'Input') && isfield(e.Input,'cluster') && isfield(e.Input.cluster,'wall_time')
        w = string(e.Input.cluster.wall_time);
    end
    e.request = struct('wall_time', w, 'num_workers', NaN, 'num_threads', NaN, 'partition', "");
end
e = orderfields(e, {'key','type','profile','job_id','delete_job','submitted','fullpath', ...
                    'filename','kind','plot_on_collect','request','Input'});
end

function t = job_time(h, names)
% First available datetime property of a job or future (NaT if none)
t = NaT;
for n = names
    try
        if isprop(h, char(n))
            v = h.(char(n));
            if isdatetime(v) && ~isempty(v) && ~isnat(v), t = v; return; end
        end
    catch
    end
end
end

function b = folder_bytes(loc, stem)
% Size of a job's files in the job storage: the folder <stem> and files <stem>.*
d = dir(fullfile(loc, [stem '.*']));
b = sum([d(~[d.isdir]).bytes]) + tree_bytes(fullfile(loc, stem));
end

function b = tree_bytes(folder)
% All files below folder (explicit recursion, works in every MATLAB version)
b = 0;
if ~isfolder(folder), return; end
d = dir(folder);
for n = 1:numel(d)
    if d(n).isdir
        if ~any(strcmp(d(n).name, {'.','..'})), b = b + tree_bytes(fullfile(folder, d(n).name)); end
    else
        b = b + d(n).bytes;
    end
end
end

function f = registry_folder()
if ispc
    base = getenv('LOCALAPPDATA');
else
    base = fullfile(getenv('HOME'), '.local', 'share');
end
f = fullfile(base, 'DDA_jobs');
if ~isfolder(f), mkdir(f); end
end

function f = entry_file(key)
f = fullfile(registry_folder(), [char(key) '.mat']);
end
