function E = dda_time_estimate(Input, workers)
% DDA_TIME_ESTIMATE  Estimated solve and run time of an eigen-solver run.
%
%   E = dda_time_estimate(Input)            workers from Input.parallel / Input.cluster
%   E = dda_time_estimate(Input, workers)
%
% Model (single-threaded workers, Roihu):
%   run   = job start-up + per-worker overhead + frequencies per worker x t(Ndof)
%   t(Ndof) per frequency: log-log interpolation between measured points
%   (dda_time_benchmarks below), growing as Ndof^3 beyond the largest point
%   (measured slope there: 2.9) and as Ndof^2 below the smallest.
%
% Fields of E:
%   ndof, nk, workers, freqs_per_worker
%   per_freq             s per frequency for one worker (eig of an Ndof x Ndof matrix)
%   per_worker_overhead  s, once per worker (warm-up, receiving the inputs)
%   job_startup          s, starting the worker pool (cluster/batch jobs)
%   solve                s, the solve as dda_solve_job times it
%   run                  s, solve + job start-up (what counts against the wall time)
%   layers               the layer benchmark (dda_layer_benchmark) with per_freq_used,
%                        the s per frequency added for building the matrix; [] without layers
%   layer_note           one line on how the layer cost was obtained ("" without layers)
%
% Used by DDA_parameter_GUI (estimates) and collect_dda_jobs (measured vs
% estimated). To add a measurement, add a row to dda_time_benchmarks: the
% per-frequency time collect_dda_jobs prints for a finished job.

switch string(Input.dimensionality)
    case "3D", nc = 3;
    case {"2D","cyl_TE"}, nc = 2;
    otherwise, nc = 1;
end
N = prod(Input.arrayN);
if isfield(Input, 'array_type') && string(Input.array_type) == "double_rectangular", N = 2*N; end
E.ndof = nc*N;
E.nk   = numel(Input.range.real_parts)*numel(Input.range.imag_parts);

if nargin < 2 || isempty(workers)
    switch string(getf(Input, 'parallel', "none"))
        case "cluster", workers = getf(getf(Input, 'cluster', struct()), 'num_workers', 1);
        case "threads", workers = getf(getf(Input, 'threads', struct()), 'num_workers', []);
                        if isempty(workers), workers = NaN; end
        otherwise,      workers = 1;
    end
end
E.workers = workers;
E.freqs_per_worker = ceil(E.nk/max(1, workers));
if ~isfinite(workers), E.freqs_per_worker = E.nk; end

B = dda_time_benchmarks();
ref = B(:,1).';  tref = B(:,2).';
n = E.ndof;
if n <= ref(1),       t = tref(1)*(n/ref(1))^2;
elseif n >= ref(end), t = tref(end)*(n/ref(end))^3;
else,                 t = exp(interp1(log(ref), log(tref), log(n)));
end
% Layered background (cyl_TM / cyl_TE with Input.layers.enable): building the
% matrix adds Sommerfeld integrals per frequency. The build is timed on this
% computer for the run's own array and layer stack (dda_layer_benchmark,
% cached) and added to the eig model. If that is switched off or fails, a
% rough model: it grows as the square of the array length (distinct
% separations x integration points), about 0.1 s (TM) and 0.15 s (TE) for 184
% cylinders on a laptop.
E.layers = [];
E.layer_note = "";
if isfield(Input, 'layers') && isfield(Input.layers, 'enable') && any(Input.layers.enable)
    Bl = dda_layer_benchmark(Input);
    if Bl.ok
        t_layers = Bl.per_freq;
    else
        if nc == 2, c = 0.15; else, c = 0.10; end
        t_layers = c*(N/184)^2;
    end
    t = t + t_layers;
    E.layers = Bl;
    E.layers.per_freq_used = t_layers;
    E.layer_note = Bl.note;
end
E.per_freq = t;
E.per_worker_overhead = 0.7;
E.job_startup = 6;
E.solve = E.per_worker_overhead + E.freqs_per_worker*t;
E.run   = E.solve + E.job_startup;
end

function B = dda_time_benchmarks()
% [Ndof, s per frequency per single-threaded worker], fixed costs removed:
% (solve - 0.7 s) / frequencies per worker. Measured on Roihu, complex
% frequencies (lean results, no tracking).
B = [  200   0.03     % 1Dx 20x10, job 2910437: 1.0 s, 10 freqs/worker (uncertain: mostly overhead)
       400   0.090    % 1Dx 20x20, job 2913916: 181 s, 2000 freqs/worker
       800   0.593    % 2D  20x20, job 2918925: 238 s, 400 freqs/worker
      1200   2.83     % 3D  20x20, early benchmark 2.9 s at ~10 freqs/worker
      1800  10.5      % 2D  30x30, job 2910231: 106 s, 10 freqs/worker
      3200  56.8 ];   % 2D  40x40, job 2919206: 1136 s, 20 freqs/worker
end

function v = getf(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
