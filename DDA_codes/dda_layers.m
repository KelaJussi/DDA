function L = dda_layers(Input, Array)
% DDA_LAYERS  The layered background of a cylinder run (cyl_TM, cyl_TE).
%
%   L = dda_layers(Input)          the stack (L.enable = false without layers)
%   L = dda_layers(Input, Array)   also the height of every cylinder, L.z (N x 1)
%
% Geometry: the cylinders lie along y, the array runs along x and the layers
% are stacked along z. The cylinders sit inside one layer, the "host layer",
% whose index is the host index of the run (Input.n or Input.substrate, as
% without layers). Any number of layers can be put below and above it:
%
%   Input.layers.enable          true to use the layered background
%   Input.layers.host_thickness  thickness of the host layer [m]
%   Input.layers.z_cyl           height of the cylinder axes above the lower
%                                interface of the host layer [m]; the z of
%                                Array.r_i is added to it (0 < z < thickness)
%   Input.layers.below_n         indices below the host layer, nearest first;
%                                the last one is the half-space (substrate)
%   Input.layers.below_d         thickness [m] of each layer below except the
%                                last: numel(below_n) - 1 values
%   Input.layers.above_n         indices above the host layer, nearest first;
%                                the last one is the half-space (cover)
%   Input.layers.above_d         thickness [m] of each layer above except the
%                                last: numel(above_n) - 1 values
%
% Example, air / 1.18 um film (host) / 2 um MY133 / glass:
%   Input.layers.host_thickness = 1.18e-6;
%   Input.layers.below_n = [1.33 1.52];   Input.layers.below_d = 2e-6;
%   Input.layers.above_n = 1;             Input.layers.above_d = [];
%
% Indices are constants and may be complex (n + i*kappa, kappa > 0 absorbs).
% A half-space with the host index removes that interface.
L = struct('enable', false);
if ~isfield(Input, 'layers') || ~isstruct(Input.layers), return; end
S = Input.layers;
if ~isfield(S, 'enable') || ~any(S.enable), return; end
if ~any(strcmp(char(Input.dimensionality), {'cyl_TM', 'cyl_TE'}))
    error('dda_layers:dimensionality', ...
        'Input.layers is supported for the cylinder models only (dimensionality "cyl_TM" or "cyl_TE").');
end

L.enable = true;
L.d  = double(need(S, 'host_thickness'));
L.nb = double(reshape(need(S, 'below_n'), 1, []));
L.na = double(reshape(need(S, 'above_n'), 1, []));
L.db = double(reshape(opt(S, 'below_d', []), 1, []));
L.da = double(reshape(opt(S, 'above_d', []), 1, []));
L.z_cyl = double(opt(S, 'z_cyl', L.d/2));

if ~isscalar(L.d) || ~(L.d > 0)
    error('dda_layers:thickness', 'Input.layers.host_thickness must be a positive length in metres.');
end
check_side(L.nb, L.db, 'below');
check_side(L.na, L.da, 'above');

if nargin > 1 && ~isempty(Array)
    r = double(Array.r_i);
    L.z = L.z_cyl + r(:,3);
    if any(L.z <= 0) || any(L.z >= L.d)
        error('dda_layers:position', ['Every cylinder axis must lie inside the host layer: ' ...
            '0 < z_cyl + z < host_thickness (z from %.1f to %.1f nm, thickness %.1f nm).'], ...
            min(L.z)*1e9, max(L.z)*1e9, L.d*1e9);
    end
    if any(abs(r(:,2) - r(1,2)) > 1e-12)
        error('dda_layers:y', ['With layers the cylinders (along y) must all have the same y: ' ...
            'the array runs along x (use arrayN = [Nx 1 1]) and the layers are stacked along z.']);
    end
end
end

function check_side(n, d, name)
if isempty(n)
    error('dda_layers:n', 'Input.layers.%s_n needs at least one index (the half-space).', name);
end
if numel(d) ~= numel(n) - 1
    error('dda_layers:d', ['Input.layers.%s_d needs one thickness for every layer except the last ' ...
        '(the half-space): %d value(s) for %d indices, got %d.'], name, numel(n)-1, numel(n), numel(d));
end
if any(d <= 0) || any(real(n) <= 0)
    error('dda_layers:values', 'Input.layers.%s: thicknesses and the real part of the indices must be positive.', name);
end
end

function v = need(S, name)
if ~isfield(S, name) || isempty(S.(name))
    error('dda_layers:missing', 'Input.layers.%s is required when Input.layers.enable is true.', name);
end
v = S.(name);
end

function v = opt(S, name, default)
if isfield(S, name) && ~isempty(S.(name)), v = S.(name); else, v = default; end
end
