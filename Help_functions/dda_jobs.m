function dda_jobs(action, sel)
% DDA_JOBS  List, cancel or remove submitted DDA jobs (without fetching).
%
%   dda_jobs                  list all jobs, numbered, with their state
%   dda_jobs("cancel", n)     stop job(s) n and forget them; n = number(s) from
%                             the list, or "all"
%   dda_jobs("remove", n)     forget job(s) n, e.g. failed ones, so the run can be
%                             submitted again (a cluster job's data is deleted)
%
% To fetch and save finished jobs, use collect_dda_jobs.

if nargin < 1, action = "list"; end
entries = dda_job_registry("list");
if isempty(entries)
    fprintf('No DDA jobs submitted.\n');
    return
end

switch lower(string(action))
    case "list"
        fprintf('%4s  %-9s %-10s %-19s %s\n', 'n', 'state', 'where', 'submitted', 'result');
        for n = 1:numel(entries)
            e = entries(n);
            try, st = dda_job_registry("state", e); catch, st = "?"; end
            where = e.type;
            if e.type == "cluster", where = sprintf('job %d', e.job_id); end
            fprintf('%4d  %-9s %-10s %-19s %s\n', n, st, where, ...
                char(e.submitted, 'yyyy-MM-dd HH:mm'), e.filename);
        end
    case {"cancel","remove"}
        if nargin < 2, error('Give the job number(s) from dda_jobs, or "all".'); end
        if isstring(sel) || ischar(sel)
            if lower(string(sel)) ~= "all", error('Use job numbers or "all".'); end
            sel = 1:numel(entries);
        end
        for n = sel(:).'
            if n < 1 || n > numel(entries)
                fprintf('No job %d.\n', n); continue
            end
            e = entries(n);
            if lower(string(action)) == "cancel"
                try, dda_job_registry("cancel", e); catch, end
            end
            dda_job_registry("remove", e);
            if lower(string(action)) == "cancel", done = "cancelled"; else, done = "removed"; end
            fprintf('Job %d %s: %s\n', n, done, e.filename);
        end
    otherwise
        error('dda_jobs:action', 'Use dda_jobs, dda_jobs("cancel", n) or dda_jobs("remove", n).');
end
end
