function filename = construct_filename(Input, kind)
% Build a descriptive .mat filename from the simulation Input struct.
%
%   construct_filename(Input)          complex-frequency results (as before)
%   construct_filename(Input, "real")  real-frequency results: adds the
%                                      material, incident field and solver,
%                                      and gives the range in nm
%
% The name is a comma-separated list of parts. Each part comes from one
% small helper function below; a helper returns '' when its part does not
% apply, and empty parts are dropped. To add a new section to the name,
% write a helper and add one line to the parts list.
%
% For configurations without the new options, the name is identical to
% the one from the previous version, so old results can still be found.

if nargin < 2, kind = "complex"; end
is_real = lower(string(kind)) == "real";

parts = {
    char(Input.dimensionality)
    sprintf('%d x %d', Input.arrayN(1), Input.arrayN(2))
    sprintf('p=(%d,%d)', floor(Input.p(1)*1e9), floor(Input.p(2)*1e9))
    part_array_type(Input)
    part_periodic(Input)      % new: only added when a periodic boundary is on
    part_HG(Input)
    sprintf('r=%d nm', ceil(Input.particle.r*1e9))
    part_particle(Input)      % new: only for non-Mie particles (Lorentzian)
    part_array_edit(Input)    % new: only when the layout was edited in the array editor
    part_cylinders(Input)     % new: only for cyl_TM / cyl_TE (materials, layered background)
    };
if is_real
    parts = [parts; {
        part_material(Input)
        part_field(Input)
        part_solver(Input)
        part_range_real(Input)}];
else
    parts = [parts; {part_range(Input)}];
end
parts = [parts; {char(Input.save_message)}];

parts = parts(~cellfun(@isempty, parts));
filename = string([strjoin(parts.', ', ') '.mat']);
end

%% ---------------------------------------------------------------- parts
function s = part_array_type(Input)
% Section that depends on the array type
switch char(Input.array_type)
    case 'double_rectangular'
        o = Input.second_array_offset*1e9;
        s = sprintf('double [%.4f, %.4f, %.4f]', o(1), o(2), o(3));
    case 'PC_mirrors'
        s = sprintf('PC %d NPs [%.3e, %.3e; %.3e, %.3e]', Input.PC_mirror_thickness, ...
            Input.PC_mirror_p(1), Input.PC_mirror_p(2), ...
            Input.PC_mirror_offset(1), Input.PC_mirror_offset(2));
    otherwise
        m = Input.mirror_A;
        s = sprintf('mirror [%d,%d,%d,%d]', m(1), m(2), m(3), m(4));
end
end

function s = part_periodic(Input)
% Periodic boundaries, e.g. 'PBC [1,0] L=(300,0) nm Nc=(10,0)'
% (same defaults as Make_int_matrix_A_SI; 'auto' = cell size inferred
% from the particle positions)
s = '';
if ~isfield(Input,'periodic_boundary') || ~any(Input.periodic_boundary)
    return
end
pbc = logical(Input.periodic_boundary(:).');

Nc = [10 10];
if isfield(Input,'periodic_Ncells')
    Nc = double(Input.periodic_Ncells(:).') .* [1 1];
end
Nc(~pbc) = 0;

L = [0 0];
if isfield(Input,'periodic_cell') && ~isempty(Input.periodic_cell)
    L = double(Input.periodic_cell(:).') .* [1 1];
end
L(~pbc) = 0;
Lstr = {sprintf('%d', round(L(1)*1e9)), sprintf('%d', round(L(2)*1e9))};
Lstr(pbc & L <= 0) = {'auto'};

s = sprintf('PBC [%d,%d] L=(%s,%s) nm Nc=(%d,%d)', pbc(1), pbc(2), ...
    Lstr{1}, Lstr{2}, Nc(1), Nc(2));
end

function s = part_HG(Input)
if strcmp(Input.array_type, "1Dsquare_halfshifts2")
    halfshifttype = 2;
else
    halfshifttype = 1;
end
s = sprintf('HG %d%d (%d)', Input.array_shifts_y, Input.array_shifts_x, halfshifttype);
end

function s = part_range(Input)
% Kept exactly as before (incl. the missing space in '(N=..)eV') so that
% names of old result files still match.
re = Input.range.real_parts;
im = Input.range.imag_parts;
s = sprintf('range real %.4f - %.4f (N=%d)eV, imag %.4f - %.4f (N=%d) eV', ...
    re(1), re(end), length(re), im(1), im(end), length(im));
end
%% ------------------------------------------- parts for real-frequency runs
function s = part_material(Input)
% Material and substrate; the radii only when the particle is anisotropic
% (Lorentzian particles are described by part_particle instead)
s = '';
if is_lorentzian(Input), return; end
s = sprintf('%s on %s', Input.sample.type, Input.substrate);
if string(Input.sample.type) == "glass_Mie"
    s = sprintf('%s n=%.3f', s, Input.Mie_glass_n);
end
if isfield(Input,'particle_gain') && Input.particle_gain ~= 0
    s = sprintf('%s gain=%g', s, Input.particle_gain);
end
if isfield(Input.particle,'anisotropic') && Input.particle.anisotropic
    s = sprintf('%s r=(%d,%d,%d) nm', s, round(Input.particle.rx*1e9), ...
        round(Input.particle.ry*1e9), round(Input.particle.rz*1e9));
end
end

function s = part_field(Input)
% Incident field: polarization, angles and the parameters of the preset
s = sprintf('Ei %s th=%g ph=%g psi=%g', Input.polarization, Input.theta, Input.phi, Input.psi);
preset = string(Input.preset_field);
switch preset
    case "none"
        return
    case {"sine_y","exp_y"}
        s = sprintf('%s %s ord=%g ph=%.3f', s, preset, Input.preset_sine_order, Input.preset_sine_phase);
    case "gaussian_proper"
        s = sprintf('%s %s w0=%.2f um z=%.2f um HG%d%d', s, preset, Input.Gaussian_w0*1e6, ...
            Input.Gaussian_z*1e6, Input.Hermite_order_x, Input.Hermite_order_y);
    case {"hermite_gaussian","hermite_JTK_tests"}
        s = sprintf('%s %s HG%d%d sc=(%g,%g) off=(%g,%g) R=%g', s, preset, ...
            Input.Hermite_order_x, Input.Hermite_order_y, Input.Hermite_scale_x, ...
            Input.Hermite_scale_y, Input.Hermite_offset_x, Input.Hermite_offset_y, Input.Hermite_R);
    otherwise
        s = sprintf('%s %s', s, preset);
end
end

function s = part_solver(Input)
% Eigen-solver files also hold the mode data the plots need
if string(Input.inversion_tool) == "eig", s = 'eig'; else, s = 'direct'; end
end

function s = part_range_real(Input)
wl = Input.range.wl*1e9;
s = sprintf('wl %.2f - %.2f nm (N=%d)', min(wl), max(wl), numel(wl));
end

%% ------------------------------------------------ edited layout (both kinds)
function s = part_array_edit(Input)
% '' without edits, so names of unedited runs do not change
s = '';
if ~isfield(Input, 'array_edit') || isempty(Input.array_edit) || ~Input.array_edit.enable, return; end
E = Input.array_edit;
s = sprintf('edited layout (%d moved, %d resized) #%s', E.n_moved, E.n_resized, E.hash);
end

%% ------------------------------------------------ cylinders and layers (both kinds)
function s = part_cylinders(Input)
% '' unless the run uses cylinders. Materials of the cylinders and the host
% (complex-frequency names do not have them otherwise) and the layer stack,
% e.g. 'cyl n=1.33 in 1.4665, layers 1+host 1180nm z100+1.33 2000nm+1.52'
% (top to bottom; no characters that Windows file names forbid)
s = '';
if ~any(strcmp(char(Input.dimensionality), {'cyl_TM','cyl_TE'})), return; end
if strcmp(char(Input.sample.type), 'glass_Mie')
    s = sprintf('cyl n=%s', num2str(Input.Mie_glass_n, 6));
    if isfield(Input,'particle_gain') && Input.particle_gain ~= 0
        s = sprintf('%s gain=%g', s, Input.particle_gain);
    end
else
    s = sprintf('cyl %s', char(Input.sample.type));
end
if strcmp(char(Input.substrate), 'const'), s = sprintf('%s in %s', s, num2str(Input.n, 6));
else,                                      s = sprintf('%s in %s', s, char(Input.substrate));
end
if isfield(Input.particle,'cyl_alpha_legacy') && Input.particle.cyl_alpha_legacy
    s = [s ' old alpha'];
end
if ~isfield(Input,'layers') || ~isfield(Input.layers,'enable') || ~Input.layers.enable, return; end
Ly = Input.layers;
zc = Ly.host_thickness/2;
if isfield(Ly,'z_cyl') && ~isempty(Ly.z_cyl), zc = Ly.z_cyl; end
above = stack_str(Ly.above_n, getfield_or(Ly,'above_d'));
below = stack_str(Ly.below_n, getfield_or(Ly,'below_d'));
s = sprintf('%s, layers %s+host %gnm z%g+%s', s, strjoin(fliplr(above), '+'), ...
    Ly.host_thickness*1e9, zc*1e9, strjoin(below, '+'));
end

function c = stack_str(n, d)
% {'1.33 2000nm', '1.52'}: nearest layer first, the half-space has no thickness
c = cell(1, numel(n));
for m = 1:numel(n)
    c{m} = num2str(n(m), 6);
    if m <= numel(d), c{m} = sprintf('%s %gnm', c{m}, d(m)*1e9); end
end
end

function v = getfield_or(S, name)
if isfield(S, name), v = S.(name); else, v = []; end
end

%% ------------------------------------------------ particle model (both kinds)
function s = part_particle(Input)
% '' for the Mie particles, so complex names of Mie runs stay as before.
% Lorentzian (MJH) particles: substrate, A0, tau, peaks, units and correction.
s = '';
if ~is_lorentzian(Input), return; end
P = Input.particle;
units = 'gaussian'; if isfield(P,'lorentz_units'), units = char(P.lorentz_units); end
corr  = 'none';     if isfield(P,'lorentz_correction'), corr = char(P.lorentz_correction); end
s = sprintf('Lorentzian on %s A0=%s tau=%s fs peak=%s nm %s %s', Input.substrate, ...
    fmt_values(P.scatter_strength, 1, '%.3g'), fmt_values(P.tau, 1e15, '%.3g'), ...
    fmt_values(P.lambda_peak_vec, 1e9, '%.1f'), units, corr);
if strcmpi(corr, 'MLWA') && isfield(P,'lorentz_mlwa_a')
    s = sprintf('%s a=%d nm', s, round(P.lorentz_mlwa_a*1e9));
end
end

function tf = is_lorentzian(Input)
tf = any(strcmpi(string(Input.sample.type), ["Lorentzian","lorentzian Au","nonidentical"]));
end

function s = fmt_values(v, scale, f)
% 'a' if all equal, '[a,b,c]' for up to 3 values, else size and a checksum
v = double(v(:).')*scale;
if all(v == v(1))
    s = sprintf(f, v(1));
elseif numel(v) <= 3
    s = ['[' strjoin(arrayfun(@(x) sprintf(f, x), v, 'UniformOutput', false), ',') ']'];
else
    h = mod(sum(round(v/max(abs(v))*1e6).*(1:numel(v))), 2^31);
    s = sprintf('%d values #%x', numel(v), h);
end
end
