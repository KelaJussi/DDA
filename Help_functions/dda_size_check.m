function dda_size_check(Input, Array)
% DDA_SIZE_CHECK  Stop before submitting a job that cannot fit, and say why.
%
% Called by submit_dda_job. Limits checked:
%   * 2 GB per variable in a .mat file: batch saves the job's inputs (Input,
%     Array) and outputs (the results) in .mat files; complex-frequency
%     results are also saved by save_complex_freqs_modes in that format.
%   * Cluster memory per CPU (Input.cluster.mem_per_cpu x num_threads): the
%     job's controlling worker holds the inputs and all results; each pool
%     worker holds a copy of the inputs and one eigenproblem at a time.
% Prints a one-line summary when everything fits.

GB   = 1e9;
lim  = 1.9e9;                                  % 2 GB per variable, with a margin
w = whos('Input');  in_b  = w.bytes;
w = whos('Array');  arr_b = w.bytes;
R = dda_result_size(Input);
mode = string(Input.parallel);

problems = strings(0);
if arr_b > lim
    problems(end+1) = sprintf(['Array (positions, polarizabilities) is %.2f GB; job inputs are saved in a ' ...
        '.mat file, limited to 2 GB per variable.'], arr_b/GB);
end
if in_b > lim
    problems(end+1) = sprintf('Input is %.2f GB, over the 2 GB per-variable limit for job inputs.', in_b/GB);
end
if R.results > lim
    problems(end+1) = sprintf(['The results would be about %.2f GB, over the 2 GB per-variable limit ' ...
        'for returning and saving them.'], R.results/GB);
end

mem = NaN;
if mode == "cluster" && isfield(Input, 'cluster') && ~strcmpi(string(getf(Input.cluster, 'name', "")), "Processes")
    mem = parse_mem(getf(Input.cluster, 'mem_per_cpu', "")) * max(1, getf(Input.cluster, 'num_threads', 1));
    if isfinite(mem)
        client = in_b + arr_b + R.results + R.vectors;
        worker = in_b + arr_b + R.eig_workspace;
        if client > mem
            problems(end+1) = sprintf(['The job''s controlling worker needs about %.1f GB (inputs, all results%s), ' ...
                'but gets %.1f GB (mem_per_cpu x num_threads).'], client/GB, ...
                ternary(R.track, ', eigenvectors for tracking', ''), mem/GB);
        end
        if worker > mem
            problems(end+1) = sprintf(['Each pool worker needs about %.1f GB (a copy of the inputs and one ' ...
                '%d x %d eigenproblem), but gets %.1f GB.'], worker/GB, R.ndof, R.ndof, mem/GB);
        end
    end
end

if isempty(problems)
    msg = sprintf('Sizes: inputs %s, results ~%s', human(in_b + arr_b), human(R.results));
    if R.lean, msg = msg + " (lean complex results)"; end
    if isfinite(mem), msg = msg + sprintf(', memory per worker %s', human(mem)); end
    fprintf('%s.\n', msg);
    return
end

% Advice: how many frequencies per run would fit
per_freq = (R.results + R.vectors)/max(R.nk, 1);      % what grows with the number of frequencies
budget = lim;
if isfinite(mem), budget = min(budget, mem - in_b - arr_b - R.eig_workspace); end
nk_max = floor(max(budget, 0)/max(per_freq, 1));
advice = "";
if nk_max >= 1 && nk_max < R.nk
    nre = numel(Input.range.real_parts);
    advice = sprintf(['\nAbout %d frequencies fit in one run: split the grid into %d runs, for example ' ...
        '%d imaginary points per run with all %d real points.'], nk_max, ceil(R.nk/nk_max), ...
        max(1, floor(nk_max/nre)), nre);
end
if R.track
    advice = advice + newline + "Tracking keeps every eigenvector: Input.track.enable = false reduces this a lot.";
end
if ~R.lean && any(Input.range.imag_parts ~= 0)
    advice = advice + newline + "Input.save.complex_all = false keeps only what the complex plots use.";
end
error('dda_size_check:tooLarge', ['This run does not fit (%d x %d matrix, %d frequencies):\n  - %s%s'], ...
    R.ndof, R.ndof, R.nk, strjoin(problems, newline + "  - "), advice);
end

%% ------------------------------------------------------------------------
function b = parse_mem(s)
% Slurm memory: "4g", "4gb", "4000m", "4096M", plain number = MB. NaN if unset.
b = NaN;
s = strtrim(lower(string(s)));
if strlength(s) == 0, return; end
t = regexp(s, '^([\d.]+)\s*([kmgt]?)b?$', 'tokens', 'once');
if isempty(t), return; end
f = struct('k', 1e3, 'm', 1e6, 'g', 1e9, 't', 1e12);
u = t{2}; if u == "", u = "m"; end
b = str2double(t{1}) * f.(char(u));
end

function v = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

function s = human(b)
if b >= 1e9, s = sprintf('%.2f GB', b/1e9); else, s = sprintf('%.0f MB', b/1e6); end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
