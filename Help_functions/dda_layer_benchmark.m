function B = dda_layer_benchmark(Input)
% DDA_LAYER_BENCHMARK  Measured cost of the layered background per frequency.
%
%   B = dda_layer_benchmark(Input)
%
% With Input.layers.enable (cyl_TM / cyl_TE), every frequency adds the
% Sommerfeld integrals of dda_layer_green to building the interaction matrix.
% Their cost depends on the array (the number of distinct separations and
% heights, and the length of the array), the layer stack and the frequency,
% so it is measured here instead of modelled: the run's own array is built,
% and the run's matrix builder (Make_cyl_TM/TE_matrix_A_SI) is timed with a
% stand-in polarizability at the corners of the frequency range (lowest and
% highest real part, at the most negative imaginary part, where the
% integration path is longest). The eig itself stays with the Roihu model of
% dda_time_estimate.
%
% Timing is single-threaded (like the cluster workers). For cluster runs the
% local time is scaled by the speed of this computer relative to Roihu,
% measured once per MATLAB session with an eig of the size of a Roihu
% benchmark row (dda_time_estimate). Results are cached per geometry, layer
% stack and frequency range, so the GUI can call this on every refresh.
%
% Fields of B:
%   ok            true if the measurement worked (otherwise see B.note)
%   per_freq      s per frequency to build the matrix, for the run's workers
%   local         the same, measured on this computer (1 thread)
%   speed_factor  Roihu time / local time (1 for local runs)
%   samples       number of frequencies timed
%   cached        true if taken from the cache
%   note          one line for the GUI
%
% Options (all optional):
%   Input.layers.benchmark         false: skip, dda_time_estimate uses a rough model
%   Input.layers.benchmark_budget  s spent timing at most (default 5); with a
%                                  slow first sample only that one is used
%   Input.layers.speed_factor      fixed Roihu/local factor instead of measuring it

persistent cache ref_speed
B = struct('ok', false, 'per_freq', NaN, 'local', NaN, 'speed_factor', 1, ...
           'samples', 0, 'cached', false, 'note', "");
if ~isfield(Input, 'layers') || ~isfield(Input.layers, 'enable') || ~any(Input.layers.enable)
    return
end
S = Input.layers;
if isfield(S, 'benchmark') && ~isempty(S.benchmark) && ~any(S.benchmark)
    B.note = "layer cost not benchmarked (layers.benchmark = false): rough model";
    return
end
budget = 5;
if isfield(S, 'benchmark_budget') && ~isempty(S.benchmark_budget), budget = S.benchmark_budget; end
on_cluster = strcmp(char(string(getf(Input, 'parallel', "none"))), 'cluster');

key = cache_key(Input);
if isempty(cache), cache = struct('key', {}, 'B', {}); end
hit = find(strcmp({cache.key}, key), 1);
if ~isempty(hit)
    B = cache(hit).B;
    B.cached = true;
    B = apply_speed(B, on_cluster, S, ref_speed);
    return
end

try
    %% The run's array (positions), as the run builds it
    try
        [Array, Ib] = Make_generic_Array(Input);
        if exist('dda_apply_array_edit', 'file')
            [Array, ~] = dda_apply_array_edit(Array, Ib);
        end
        r = double(Array.r_i);
    catch
        nx = Input.arrayN(1);                            % fallback: a straight chain along x
        r  = [(0:nx-1).'*Input.p(1), zeros(nx, 2)];
        Array = struct('r_i', r);
    end
    L = dda_layers(Input, Array);                        % checks the stack and positions
    N = size(r, 1);
    te = strcmp(char(string(Input.dimensionality)), 'cyl_TE');
    nc = 1 + te;
    nh = host_index(Input);

    %% Frequencies: the corners where the integral is most expensive
    hc  = 1239.84187e-9;                                 % eV*m
    re  = Input.range.real_parts;  im = Input.range.imag_parts;
    eVs = unique([max(re), min(re)] + 1i*min(im), 'stable');
    k0s = 2*pi*eVs/hc;

    %% The run's matrix builder with a stand-in polarizability (its value does
    %  not change the cost), so the timing includes everything the build does
    Inb = Input;  Inb.N = N;
    Arb = struct('r_i', r, 'k', nh*k0s, 'n', nh*ones(size(k0s)), 'alpha_type', ones(N, 1));
    Arb.alpha_types = zeros(1, nc, nc, numel(k0s));
    for c = 1:nc, Arb.alpha_types(1, c, c, :) = 1e-30*(1 + 0.1i); end
    if te, build = @(kk) Make_cyl_TE_matrix_A_SI(Inb, Arb, kk);
    else,  build = @(kk) Make_cyl_TM_matrix_A_SI(Inb, Arb, kk);
    end

    %% Time it, single-threaded
    restore = single_thread(); %#ok<NASGU>  restores the thread count when cleared
    if te, pol = 'TE'; else, pol = 'TM'; end
    dda_layer_green(k0s(1), nh, L, [0 1e-6], L.z([1 1]).', L.z([1 1]).', pol);   % warm-up
    t = zeros(1, 0);
    spent = tic;
    for m = 1:numel(k0s)
        t0 = tic;
        A = build(m); %#ok<NASGU>
        t(end+1) = toc(t0); %#ok<AGROW>
        if toc(spent) > budget, break; end
    end
    clear A
    B.local   = mean(t);
    B.samples = numel(t);
    B.ok      = true;
    B.note    = sprintf("matrix build with layers benchmarked: %.3g s per frequency on this computer (%d cylinders, %d sample%s)", ...
                        B.local, N, B.samples, plural(B.samples));
catch err
    B.note = "layer benchmark failed (" + string(err.message) + "): rough model";
    cache(end+1) = struct('key', key, 'B', B);      % do not retry on every GUI refresh
    return
end

%% Speed of this computer relative to Roihu, once per session
if on_cluster && isempty(ref_speed) && ~(isfield(S, 'speed_factor') && ~isempty(S.speed_factor))
    ref_speed = measure_speed();
end
cache(end+1) = struct('key', key, 'B', B);
if numel(cache) > 50, cache = cache(end-49:end); end
B = apply_speed(B, on_cluster, S, ref_speed);
end

%% ------------------------------------------------------------------------
function B = apply_speed(B, on_cluster, S, ref_speed)
if ~B.ok, return; end
f = 1;
if on_cluster
    if isfield(S, 'speed_factor') && ~isempty(S.speed_factor), f = S.speed_factor;
    elseif ~isempty(ref_speed) && isfinite(ref_speed),         f = ref_speed;
    end
end
B.speed_factor = f;
B.per_freq = B.local*f;
if on_cluster && f ~= 1
    B.note = B.note + sprintf(", x%.2g for Roihu", f);
end
end

function f = measure_speed()
% Roihu time / local time for eig of a 800 x 800 complex matrix, single-threaded
% (the 800 row of dda_time_estimate: 0.593 s per frequency, nearly all eig).
% NaN if it cannot be measured (then no scaling).
f = NaN;
try
restore = single_thread(); %#ok<NASGU>
rs = RandStream('mt19937ar', 'Seed', 1);
A = complex(rand(rs, 800), rand(rs, 800));
eig(A(1:50,1:50));                                      % warm-up
t = inf;
for it = 1:2
    t0 = tic;  eig(A);  t = min(t, toc(t0));
end
f = min(max(0.593/t, 0.25), 4);                         % keep it sane
catch
end
end

function c = single_thread()
% Single-threaded until the returned object is cleared (MATLAB only)
c = [];
if exist('maxNumCompThreads', 'builtin') || exist('maxNumCompThreads', 'file')
    old = maxNumCompThreads(1);
    c = onCleanup(@() maxNumCompThreads(old));
end
end

function n = host_index(Input)
% Host index of the layer with the cylinders (only its size matters here)
n = 1.5;
if strcmp(char(string(getf(Input, 'substrate', 'const'))), 'const') && isfield(Input, 'n')
    n = double(Input.n(1));
end
end

function key = cache_key(Input)
% Everything the cost depends on
f = {'dimensionality', 'arrayN', 'p', 'array_type', 'array_shifts_x', 'array_shifts_y', ...
     'array_edit', 'mirror_A', 'layers', 'n', 'substrate', 'parallel'};
parts = cell(1, numel(f));
for m = 1:numel(f)
    if isfield(Input, f{m}), parts{m} = [f{m} '=' to_str(Input.(f{m}))]; end
end
re = Input.range.real_parts;  im = Input.range.imag_parts;
parts{end+1} = sprintf('range=%.3g,%.3g,%.3g', min(re), max(re), min(im));
key = strjoin(parts(~cellfun('isempty', parts)), ';');
end

function s = to_str(v)
if isstruct(v)
    fn = fieldnames(v);  c = cell(1, numel(fn)*numel(v));  q = 0;
    for e = 1:numel(v)
        for m = 1:numel(fn)
            if any(strcmp(fn{m}, {'benchmark', 'benchmark_budget', 'speed_factor'})), continue; end
            q = q + 1;  c{q} = [fn{m} ':' to_str(v(e).(fn{m}))];
        end
    end
    s = ['{' strjoin(c(1:q), ',') '}'];
elseif iscell(v)
    c = cellfun(@to_str, v, 'UniformOutput', false);
    s = ['{' strjoin(c, ',') '}'];
elseif isnumeric(v) || islogical(v)
    if numel(v) > 64, s = sprintf('[%d:%.10g:%.10g]', numel(v), sum(double(v(:))), sum(double(v(:)).^2));
    else,             s = mat2str(double(v), 10);
    end
elseif ischar(v)
    s = v;
else
    s = strjoin(cellstr(string(v)), ',');
end
end

function v = getf(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

function s = plural(n)
if n == 1, s = ''; else, s = 's'; end
end
