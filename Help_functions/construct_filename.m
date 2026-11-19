function filename = construct_filename(Input, kind)
% Build a compact, descriptive .mat filename from the simulation Input struct.
%
%   construct_filename(Input)          complex-frequency results
%   construct_filename(Input, "real")  real-frequency results: adds the
%                                      material, incident field and solver,
%                                      and gives the range as wavelengths
%
% Fields are separated by single spaces and contain no spaces themselves.
% Units are implied: lengths in nm, except the beam waist and focus distance of
% a Gaussian field (marked um); energies in eV. Example (real frequencies):
%   1Dx 1x120 p(700,700) m[0000] PBCx700n10 HG10 r65 AuDM-BK7 L[0,0,90]Gau[80um,HG01] eig wl688.80-1239.84(N1500).mat
%   dimensionality, array size, periods, mirrors, periodic boundary (x: cell 700,
%   10 image cells per side), half-period shifts (HGyx), radius, particle on
%   substrate, field (polarization L/RC/LC [theta,phi,psi], Gaussian [w0,(z,)HGxy]),
%   solver, wavelength range (points).
% Each field comes from one small helper below; a helper returns '' when its
% field does not apply, and empty fields are dropped. Names that would make the
% path too long for Windows are shortened (fit_path).

if nargin < 2, kind = "complex"; end
is_real = lower(string(kind)) == "real";

parts = {
    char(Input.dimensionality)
    sprintf('%dx%d', Input.arrayN(1), Input.arrayN(min(2,end)))
    sprintf('p(%d,%d)', floor(Input.p(1)*1e9), floor(Input.p(2)*1e9))
    part_array_type(Input)
    part_periodic(Input)      % only when a periodic boundary is on
    part_HG(Input)            % only with half-period shifts
    sprintf('r%d', ceil(Input.particle.r*1e9))
    part_particle(Input)      % only for Lorentzian particles
    part_array_edit(Input)    % only when the layout was edited in the array editor
    part_cylinders(Input)     % only for cyl_TM / cyl_TE (materials, layered background)
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
if isfield(Input, 'save_message'), parts = [parts; {char(Input.save_message)}]; end

parts = parts(~cellfun(@isempty, parts));
filename = string([strjoin(parts.', ' ') '.mat']);
filename = fit_path(filename, kind);
end

%% ------------------------------------------------------------ path length
function filename = fit_path(filename, kind)
% Windows paths are limited to 259 characters, and MATLAB's -v7.3 (HDF5) save
% then fails with "Unable to write to file ... because it appears to be
% corrupt". If folder + name is longer, the name keeps its beginning and gets
% a fingerprint of the full name instead of the rest: "<start> ... #<hash>.mat".
% The shortening is deterministic, so a re-run finds the file again; the full
% inputs are inside the file (Input), which find_results reads.
maxlen = 255;                                    % a little below 259
folder = dda_results_folder(kind);
room   = maxlen - strlength(string(folder)) - 1;
if strlength(filename) <= room, return; end
full = char(filename);
tag  = [' ... #' dda_checksum(double(full)) '.mat'];
keep = room - numel(tag);
if keep < 20
    error('construct_filename:path', ['The results folder path is too long (%d characters) ' ...
          'for result files on Windows. Move the DDA folder to a shorter path.'], strlength(string(folder)));
end
filename = string([strtrim(full(1:keep)) tag]);
fprintf('Result name shortened to stay below the Windows path limit:\n  %s\n', filename);
end

%% ---------------------------------------------------------------- parts
function s = part_array_type(Input)
switch char(Input.array_type)
    case 'double_rectangular'
        o = Input.second_array_offset*1e9;
        s = sprintf('dbl(%g,%g,%g)', round(o(1),1), round(o(2),1), round(o(3),1));
    case 'PC_mirrors'
        s = sprintf('PC%d(%g,%g;%g,%g)', Input.PC_mirror_thickness, ...
            round(Input.PC_mirror_p(1)*1e9,1), round(Input.PC_mirror_p(2)*1e9,1), ...
            round(Input.PC_mirror_offset(1)*1e9,1), round(Input.PC_mirror_offset(2)*1e9,1));
    otherwise
        m = zeros(1,4); m(1:numel(Input.mirror_A)) = Input.mirror_A;
        if all(ismember(m, [0 1])), s = sprintf('m[%d%d%d%d]', m);       % m[1100]
        else,                       s = sprintf('m[%d,%d,%d,%d]', m);    % m[1,-1,0,0]
        end
end
end

function s = part_periodic(Input)
% e.g. PBCx700n10: periodic in x, cell 700 nm, 10 image cells on each side
% ('auto' = cell size inferred from the positions; same defaults as Make_int_matrix_A_SI)
s = '';
if ~isfield(Input,'periodic_boundary') || ~any(Input.periodic_boundary), return; end
pbc = logical(Input.periodic_boundary(:).');
Nc = [10 10];
if isfield(Input,'periodic_Ncells'), Nc = double(Input.periodic_Ncells(:).') .* [1 1]; end
L = [0 0];
if isfield(Input,'periodic_cell') && ~isempty(Input.periodic_cell)
    L = double(Input.periodic_cell(:).') .* [1 1];
end
s = 'PBC';  dirs = 'xy';
for d = find(pbc)
    if L(d) > 0, Ls = sprintf('%d', round(L(d)*1e9)); else, Ls = 'auto'; end
    s = sprintf('%s%s%sn%d', s, dirs(d), Ls, Nc(d));
end
end

function s = part_HG(Input)
% half-period shifts of the lattice (array_shifts_y, array_shifts_x); '' without
s = '';
if Input.array_shifts_y == 0 && Input.array_shifts_x == 0, return; end
s = sprintf('HG%d%d', Input.array_shifts_y, Input.array_shifts_x);
if strcmp(char(Input.array_type), '1Dsquare_halfshifts2'), s = [s 'b']; end
end

function s = part_range(Input)
re = Input.range.real_parts;
im = Input.range.imag_parts;
s = sprintf('Re%.4f-%.4f(N%d) Im%.4f-%.4f(N%d)', re(1), re(end), numel(re), im(1), im(end), numel(im));
end
%% ------------------------------------------- parts for real-frequency runs
function s = part_material(Input)
% particle-substrate, e.g. AuDM-BK7, glM1.476-n1.4665 (Lorentzian: part_particle)
s = '';
if is_lorentzian(Input), return; end
s = sprintf('%s-%s', material_abbrev(Input), substrate_str(Input));
if isfield(Input,'particle_gain') && Input.particle_gain ~= 0
    s = sprintf('%s_g%g', s, Input.particle_gain);
end
if isfield(Input.particle,'anisotropic') && Input.particle.anisotropic
    s = sprintf('%s_r(%d,%d,%d)', s, round(Input.particle.rx*1e9), ...
        round(Input.particle.ry*1e9), round(Input.particle.rz*1e9));
end
end

function a = material_abbrev(Input)
switch char(Input.sample.type)
    case 'Au_Drude_Mie',  a = 'AuDM';
    case 'Au_Mie',        a = 'AuM';
    case 'Au_Drude_MLWA', a = 'AuDMLWA';
    case 'glass_Mie',     a = sprintf('glM%s', num2str(Input.Mie_glass_n, 6));
    otherwise,            a = strrep(char(Input.sample.type), ' ', '');
end
end

function s = substrate_str(Input)
if strcmp(char(Input.substrate), 'const'), s = sprintf('n%s', num2str(Input.n, 6));
else,                                      s = char(Input.substrate);
end
end

function s = part_field(Input)
% polarization [theta,phi,psi], then the profile: Gau[w0um,(z..um,)HGxy]
switch char(Input.polarization)
    case 'linear', pol = 'L';
    otherwise,     pol = char(Input.polarization);       % RC, LC
end
s = sprintf('%s[%g,%g,%g]', pol, Input.theta, Input.phi, Input.psi);
if any(strcmpi(char(string(Input.preset_field)), {'Gaussian', 'gaussian_proper'}))
    g = sprintf('%.4gum', Input.Gaussian_w0*1e6);
    if Input.Gaussian_z ~= 0, g = sprintf('%s,z%.4gum', g, Input.Gaussian_z*1e6); end
    if Input.Hermite_order_x ~= 0 || Input.Hermite_order_y ~= 0
        g = sprintf('%s,HG%d%d', g, Input.Hermite_order_x, Input.Hermite_order_y);
    end
    s = sprintf('%sGau[%s]', s, g);
end
end

function s = part_solver(Input)
% Eigen-solver files also hold the mode data the plots need
if string(Input.inversion_tool) == "eig", s = 'eig'; else, s = 'direct'; end
end

function s = part_range_real(Input)
wl = Input.range.wl*1e9;
s = sprintf('wl%.2f-%.2f(N%d)', min(wl), max(wl), numel(wl));
end

%% ------------------------------------------------ edited layout (both kinds)
function s = part_array_edit(Input)
s = '';
if ~isfield(Input, 'array_edit') || isempty(Input.array_edit) || ~Input.array_edit.enable, return; end
E = Input.array_edit;
s = sprintf('edit(%dm,%dr)#%s', E.n_moved, E.n_resized, E.hash);
end

%% ------------------------------------------------ cylinders and layers (both kinds)
function s = part_cylinders(Input)
% '' unless the run uses cylinders: cylinder and host materials and the layer
% stack, top to bottom, e.g. cyl1.33in1.4665_lay(1+h1180z100+1.33x2000+1.52)
s = '';
if ~any(strcmp(char(Input.dimensionality), {'cyl_TM','cyl_TE'})), return; end
if strcmp(char(Input.sample.type), 'glass_Mie')
    s = sprintf('cyl%s', num2str(Input.Mie_glass_n, 6));
    if isfield(Input,'particle_gain') && Input.particle_gain ~= 0
        s = sprintf('%s_g%g', s, Input.particle_gain);
    end
else
    s = sprintf('cyl%s', material_abbrev(Input));
end
if strcmp(char(Input.substrate), 'const'), s = sprintf('%sin%s', s, num2str(Input.n, 6));
else,                                      s = sprintf('%sin%s', s, char(Input.substrate));
end
if isfield(Input.particle,'cyl_alpha_legacy') && Input.particle.cyl_alpha_legacy
    s = [s '_oldalpha'];
end
if ~isfield(Input,'layers') || ~isfield(Input.layers,'enable') || ~Input.layers.enable, return; end
Ly = Input.layers;
zc = Ly.host_thickness/2;
if isfield(Ly,'z_cyl') && ~isempty(Ly.z_cyl), zc = Ly.z_cyl; end
above = stack_str(Ly.above_n, getfield_or(Ly,'above_d'));
below = stack_str(Ly.below_n, getfield_or(Ly,'below_d'));
s = sprintf('%s_lay(%s+h%gz%g+%s)', s, strjoin(fliplr(above), '+'), ...
    Ly.host_thickness*1e9, zc*1e9, strjoin(below, '+'));
end

function c = stack_str(n, d)
% {'1.33x2000', '1.52'}: nearest layer first, the half-space has no thickness
c = cell(1, numel(n));
for m = 1:numel(n)
    c{m} = num2str(n(m), 6);
    if m <= numel(d), c{m} = sprintf('%sx%g', c{m}, d(m)*1e9); end
end
end

function v = getfield_or(S, name)
if isfield(S, name), v = S.(name); else, v = []; end
end

%% ------------------------------------------------ particle model (both kinds)
function s = part_particle(Input)
% '' for the Mie particles. Lorentzian (MJH) particles: substrate, A0, tau (fs),
% peaks, units and correction, e.g. Lor-BK7_A0=1e-07_tau=[1,4,4]fs_pk=[780,330,1000]_gaussian_none
s = '';
if ~is_lorentzian(Input), return; end
P = Input.particle;
units = 'gaussian'; if isfield(P,'lorentz_units'), units = char(P.lorentz_units); end
corr  = 'none';     if isfield(P,'lorentz_correction'), corr = char(P.lorentz_correction); end
s = sprintf('Lor-%s_A0=%s_tau=%sfs_pk=%s_%s_%s', substrate_str(Input), ...
    fmt_values(P.scatter_strength, 1, '%.3g'), fmt_values(P.tau, 1e15, '%.3g'), ...
    fmt_values(P.lambda_peak_vec, 1e9, '%g'), units, corr);
if strcmpi(corr, 'MLWA') && isfield(P,'lorentz_mlwa_a')
    s = sprintf('%s_a=%d', s, round(P.lorentz_mlwa_a*1e9));
end
end

function tf = is_lorentzian(Input)
tf = any(strcmpi(string(Input.sample.type), ["Lorentzian","lorentzian Au","nonidentical"]));
end

function s = fmt_values(v, scale, f)
% 'a' if all equal, '[a,b,c]' for up to 3 values, else count and a checksum
v = double(v(:).')*scale;
if all(v == v(1))
    s = sprintf(f, v(1));
elseif numel(v) <= 3
    s = ['[' strjoin(arrayfun(@(x) sprintf(f, x), v, 'UniformOutput', false), ',') ']'];
else
    h = mod(sum(round(v/max(abs(v))*1e6).*(1:numel(v))), 2^31);
    s = sprintf('%dvals#%x', numel(v), h);
end
end
