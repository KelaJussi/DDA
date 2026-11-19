function collected = collect_dda_jobs()
% COLLECT_DDA_JOBS  Fetch and save every DDA job that has finished.
%
%   collect_dda_jobs
%   files = collect_dda_jobs     also returns the saved result files
%
% Goes through all submitted jobs (cluster and background). Finished jobs
% are fetched (size and time are shown: the data in memory, then the saved
% file on disk), saved exactly as a blocking run would save them, removed from
% the cluster and from the job list. Queued and running jobs are listed.
% Failed jobs stay listed with their error until you remove them with
% dda_jobs("remove", n) or submit the same run again.
% Plots are drawn only for jobs submitted with Input.jobs.plot_on_collect;
% otherwise plot a saved file with real_freqs_plots_from_file(file) or
% complex_freqs_plots_from_file(name).

collected = strings(0,1);
entries = dda_job_registry("list");
if isempty(entries)
    fprintf('No DDA jobs submitted.\n');
    return
end

n_run = 0; n_queue = 0; n_fail = 0;
total_bytes = 0; total_fetch = 0;
fprintf('Checking %d DDA job(s)...\n', numel(entries));
for n = 1:numel(entries)
    e = entries(n);
    name = short_name(e.filename);
    try
        t0 = tic;
        state = dda_job_registry("state", e);   % for a remote cluster this also downloads
        t_state = toc(t0);                      % the files of a job that has just finished
    catch err
        fprintf('  [%d] %-8s could not be checked (%s): %s\n', n, "?", err.message, name);
        continue
    end
    switch state
        case "finished"
            fprintf('  [%d] finished %s\n      fetching... ', n, name);
            disk = dda_job_registry("local_bytes", e);   % job files downloaded to this computer
            try
                t0 = tic;
                R = dda_job_registry("fetch", e);
                t_fetch = toc(t0);
                times = dda_job_registry("times", e);
                R = add_timing(R, e, times);
                w = whos('R');
                t_total = t_state + t_fetch;
                if e.type == "cluster" && isfinite(disk) && disk > 0
                    fprintf(['downloaded %s in %s (%.0f Mbit/s, during the state check), ' ...
                             'read %s into memory in %s.\n'], human_size(disk), human_time(t_state), ...
                             8*disk/1e6/max(t_state, eps), human_size(w.bytes), human_time(t_fetch));
                    if t_state < 1 && 8*disk/1e6/max(t_state, eps) > 1000
                        fprintf(['      (the files were already here: an earlier dda_jobs or DDA_Master ' ...
                                 'call saw the job finish and downloaded them then)\n']);
                    end
                else
                    fprintf('%s in %s.\n', human_size(w.bytes), human_time(t_total));
                end
                R_info = [];
                if isfield(R, 'info'), R_info = R.info; end
                t0 = tic;
                finalize_dda_result(R, e, e.plot_on_collect);
                t_save = toc(t0);
                clear R
                d = dir(char(e.fullpath));
                if isempty(d), on_disk = "?"; else, on_disk = human_size(d.bytes); end
                fprintf('      saved: file %s in %s.\n', on_disk, human_time(t_save));
                print_run(e, times, R_info);
                dda_job_registry("remove", e);
                collected(end+1,1) = e.fullpath; %#ok<AGROW>
                if isfinite(disk) && disk > 0, total_bytes = total_bytes + disk;
                else,                          total_bytes = total_bytes + w.bytes;
                end
                total_fetch = total_fetch + t_total;
            catch err
                fprintf('\n      could not be saved: %s\n      (the job is kept; try again)\n', err.message);
            end
        case "failed"
            n_fail = n_fail + 1;
            [msg, diary_tail] = dda_job_registry("message", e);
            fprintf('  [%d] FAILED   %s\n      %s\n', n, name, msg);
            print_run(e, dda_job_registry("times", e), []);
            if strlength(diary_tail) > 0
                fprintf('      last output:\n%s\n', indent(diary_tail));
            end
        case "timeout"
            n_fail = n_fail + 1;
            t = dda_job_registry("times", e);
            fprintf(['  [%d] TIMEOUT  %s\n      ran past its wall time (%s): Slurm stopped the job. ' ...
                     'Raise Input.cluster.wall_time and run DDA_Master again\n      (the job is ' ...
                     'replaced automatically), or remove it with dda_jobs("remove", %d).\n'], ...
                     n, name, human_time(dda_walltime_seconds(e.request.wall_time)), n);
            if ~isnat(t.started)
                so_far = datetime('now', 'TimeZone', t.started.TimeZone) - t.started;
                fprintf('      started %s ago\n', human_time(seconds(so_far)));
            end
        case "lost"
            fprintf(['  [%d] lost     %s\n      background job of an earlier MATLAB session or ' ...
                     'of a deleted pool; removed from the list.\n'], n, name);
            dda_job_registry("remove", e);
        case "missing"
            fprintf('  [%d] missing  %s\n      the cluster job was deleted; removed from the list.\n', n, name);
            dda_job_registry("remove", e);
        case "running"
            n_run = n_run + 1;
            fprintf('  [%d] running  %s (%s)\n', n, name, age(e.submitted));
            t = dda_job_registry("times", e);
            if ~isnat(t.started)
                so_far = datetime('now', 'TimeZone', t.started.TimeZone) - t.started;
                req = parse_walltime(e.request.wall_time);
                if isnan(req)
                    fprintf('      running for %s\n', human_time(seconds(so_far)));
                else
                    fprintf('      running for %s of %s requested (%.0f%%)\n', human_time(seconds(so_far)), ...
                        human_time(seconds(req)), 100*so_far/req);
                end
            end
        otherwise
            n_queue = n_queue + 1;
            fprintf('  [%d] queued   %s (%s)\n', n, name, age(e.submitted));
    end
end
fprintf('Collected %d, running %d, queued %d, failed %d.\n', numel(collected), n_run, n_queue, n_fail);
if ~isempty(collected)
    fprintf('Fetched %s in %s in total (download and reading).\n', human_size(total_bytes), human_time(total_fetch));
end
for f = collected.'
    fprintf('  saved: %s\n', f);
end
if nargout == 0, clear collected; end          % no "ans = ..." when called as a command
end

%% ------------------------------------------------------------------------
function s = short_name(f)
s = char(f);
if numel(s) > 90, s = [s(1:60) ' ... ' s(end-24:end)]; end
end

function s = age(t)
d = datetime('now') - t;
if d < minutes(90), s = sprintf('submitted %.0f min ago', minutes(d));
elseif d < hours(48), s = sprintf('submitted %.1f h ago', hours(d));
else, s = sprintf('submitted %.0f days ago', days(d));
end
end

function R = add_timing(R, e, t)
% Scheduler times, kept with the results in C.run_info
if ~isfield(R, 'info'), R.info = struct(); end
R.info.submitted   = t.submitted;
R.info.started     = t.started;
R.info.finished    = t.finished;
R.info.queued_s    = seconds(t.started - t.submitted);
R.info.ran_s       = seconds(t.finished - t.started);
R.info.requested_wall_time = e.request.wall_time;
R.info.requested_workers   = e.request.num_workers;
R.info.requested_threads   = e.request.num_threads;
R.info.job_type    = e.type;
R.info.job_id      = e.job_id;
end

function print_run(e, t, info)
% Queue time, actual run time vs requested wall time, and what the run used
ran = t.finished - t.started;
parts = strings(0);
if ~isnat(t.started) && ~isnat(t.submitted)
    parts(end+1) = "queued " + human_time(seconds(t.started - t.submitted));
end
if ~isnan(seconds(ran))
    req = parse_walltime(e.request.wall_time);
    if isnan(req)
        parts(end+1) = "ran " + human_time(seconds(ran));
    else
        frac = ran/req;
        txt = sprintf('ran %s of %s requested (%.0f%%)', human_time(seconds(ran)), human_time(seconds(req)), 100*frac);
        if frac > 0.98,     txt = txt + " - hit the wall time limit?";
        elseif frac > 0.8,  txt = txt + " - close to the limit";
        end
        parts(end+1) = txt;
    end
end
if ~isempty(info) && isfield(info, 'solve_time_s')
    parts(end+1) = "solve " + human_time(info.solve_time_s);
end
if ~isempty(parts), fprintf('      %s\n', strjoin(parts, ', ')); end

if isempty(info), return; end                        % failed jobs: no run info
used = sprintf('%d worker(s) x %s thread(s), %s', info.workers, num_or_q(info.threads), info.device);
if info.gpus_available > 0
    used = sprintf('%s (%d GPU(s) on the node, unused)', used, info.gpus_available);
end
if e.type == "cluster" && isfinite(e.request.num_workers) && ...
        (info.workers ~= e.request.num_workers || info.threads ~= e.request.num_threads)
    used = sprintf('%s; requested %d x %s', used, e.request.num_workers, num_or_q(e.request.num_threads));
end
where = "";
if strlength(info.host) > 0,         where = " on " + info.host; end
if strlength(info.slurm_job_id) > 0, where = where + ", Slurm job " + info.slurm_job_id; end
fprintf('      %s%s\n', used, where);
print_speed(e, info);
end

function print_speed(e, info)
% Measured time per frequency per worker against the model (dda_time_estimate)
if ~isfield(info, 'solve_time_s') || ~isfield(e, 'Input') || ~isfield(e.Input, 'range'), return; end
try
    E = dda_time_estimate(e.Input, info.workers);
catch
    return
end
meas = max(info.solve_time_s - E.per_worker_overhead, 0) / E.freqs_per_worker;
fprintf('      %.3g s per frequency per worker at Ndof %d (%d frequencies per worker); model %.3g s (%+.0f%%)\n', ...
    meas, E.ndof, E.freqs_per_worker, E.per_freq, 100*(E.per_freq - meas)/max(meas, eps));
notes = strings(0);
if string(e.Input.inversion_tool) ~= "eig"
    notes(end+1) = "direct solver: the model is for the eigen solver";
elseif all(e.Input.range.imag_parts == 0)
    notes(end+1) = "real frequencies: also includes tracking, P and extinction";
end
if isfinite(info.threads) && info.threads > 1
    notes(end+1) = sprintf('%d threads per worker: the model is for 1 thread', info.threads);
end
if isfield(E, 'layers') && ~isempty(E.layers)
    notes(end+1) = sprintf('layered background: %.3g s of the model per frequency is the matrix build (%s)', ...
        E.layers.per_freq_used, E.layer_note);
end
if E.freqs_per_worker < 50
    notes(end+1) = "few frequencies per worker: the fixed ~0.7 s per worker makes this less certain";
end
if isempty(notes)
    fprintf('      benchmark row for dda_time_estimate: [%d %.3g]\n', E.ndof, meas);
else
    fprintf('      (%s)\n', strjoin(notes, '; '));
end
end

function s = num_or_q(x)
if isnan(x), s = '?'; else, s = sprintf('%d', x); end
end

function d = parse_walltime(w)
% Slurm wall time as a duration, in whole minutes as Slurm applies it
d = seconds(dda_walltime_seconds(w));
end

function s = human_size(bytes)
% Decimal units, as in file browsers: 1 MB = 1e6 bytes
if bytes >= 1e9,     s = sprintf('%.2f GB', bytes/1e9);
elseif bytes >= 1e6, s = sprintf('%.1f MB', bytes/1e6);
else,                s = sprintf('%.0f kB', bytes/1e3);
end
end

function s = human_time(sec)
if sec < 59.95, s = sprintf('%.1f s', sec);
else,           r = round(sec); s = sprintf('%d min %02d s', floor(r/60), mod(r,60));
end
end

function s = indent(t)
s = "        " + replace(t, newline, newline + "        ");
end
