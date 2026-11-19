function R = dda_solve_job(Input, Array, inner_parallel, in_cluster_job)
% DDA_SOLVE_JOB  The whole solve of one DDA run, as one function.
%
% Runs wherever the solve happens: in this MATLAB (blocking modes), on a
% background worker (parfeval) or in a cluster batch job. It returns all
% that finalize_dda_result needs to save the result file, so nothing but
% file writing is left for the client.
%
%   inner_parallel  how the eigen solver loops over frequencies:
%                   "none"    serial (background workers; parfor would be serial anyway)
%                   "threads" use the open pool, or open a thread pool (cluster
%                             batch jobs have their own pool; locally "threads")
%   in_cluster_job  true when called by a cluster batch job (submit_dda_job): the
%                   job's own pool is then never deleted. Default false.
%
%   R.P, R.C, R.Field, R.Array   solver outputs (as before)
%   R.info                       what the run used: workers, threads per worker,
%                                device, node, Slurm job id, solve time

if nargin < 4, in_cluster_job = false; end
t_start = tic;
Input.parallel = inner_parallel;
if in_cluster_job
    Input.threads.keep_pool = true;   % never delete the pool a cluster job runs in
end
% (No getCurrentJob here: it is not supported on thread-based workers.)

if Input.inversion_tool == "eig"
    if isfield(Input,'use_legacy_eig_solvers') && Input.use_legacy_eig_solvers
        error('dda_solve_job:legacy', ...
            'The legacy eigen solvers (validation) run only in this MATLAB: use Input.parallel = "none" or "threads".');
    end
    [P, C, Field, Array_out] = Calculate_DDA_eigenmodes(Input, Array);
else
    % Calculate_3D_DDA_generic opens a thread pool when useparallel is set, which
    % fails inside a worker and whenever a pool is already open; its frequency
    % loop is a plain for loop, so the pool is not used anyway.
    Input.useparallel = 0;
    Array = dda_expand_alpha(Array);   % the old builders read the full per-particle Array.alpha
    [P, C, Field, Array_out] = Calculate_3D_DDA_generic(Input, Array);
end

R = struct('P', P, 'C', C, 'Field', Field, 'Array', Array_out);
R.info = run_info(inner_parallel, toc(t_start));
end

%% ------------------------------------------------------------------------
function info = run_info(inner_parallel, t_solve)
% What this run actually used. Every query is guarded: some are not available
% on thread-based workers or without the Parallel Computing Toolbox.
info.solve_time_s   = t_solve;
info.host           = first_set(["SLURMD_NODENAME","HOSTNAME","COMPUTERNAME"]);
info.slurm_job_id   = string(getenv('SLURM_JOB_ID'));
info.workers        = 1;
info.threads        = NaN;
try, info.threads = maxNumCompThreads; catch, end
if inner_parallel == "threads"
    try
        p = gcp('nocreate');
        if ~isempty(p)
            info.workers = p.NumWorkers;
            if isa(p, 'parallel.ThreadPool'), info.threads = 1; end   % one thread per thread worker
        end
    catch
    end
end
info.device         = "CPU";      % the solvers do not use GPUs (yet)
info.gpus_available = 0;
try, if exist('gpuDeviceCount') > 0, info.gpus_available = gpuDeviceCount("available"); end, catch, end %#ok<EXIST>
info.matlab         = string(version);
end

function s = first_set(names)
s = "";
for n = names
    v = string(getenv(n));
    if strlength(v) > 0, s = v; return; end
end
end
