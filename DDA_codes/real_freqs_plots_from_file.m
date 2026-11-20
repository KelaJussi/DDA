function figs = real_freqs_plots_from_file(file, plotting)
% REAL_FREQS_PLOTS_FROM_FILE  Plots of a saved real-frequency run (replaces JTK_plots).
%
%   real_freqs_plots_from_file(filename)            options saved in Input.plotting
%   real_freqs_plots_from_file(filename, plotting)  override some of them, e.g.
%       real_freqs_plots_from_file(f, struct('emitted', struct('enable',true,'wls',1060)))
%   figs = real_freqs_plots_from_file(...)          figure handles, by plot name
%
% filename may be a full path or a name inside dda_results_folder("real"),
% or the struct returned by load(filename).
%
% Plots (Input.plotting.<name>.enable):
%   dipoles       transmission, |p| vs position and wavelength, profiles, maps;
%                 the |p| vs position panel has its own options (dipoles.position_map:
%                 energy axis, particle numbers, extinction overlay, colormap, ...)
%                 once for the driven response and once per selected eigenmode.
%                 Line arrays (1 x Ny, Nx x 1) and cylinders (cyl_TM, cyl_TE):
%                 transmission and |p| vs wavelength and position along the line
%   (Transmission is exp(-ext). If ext < 0 somewhere, the console says why:
%   usually a polarizability without radiation reaction.)
%   Q             Q-factor fit of the extinction (calc_Q_varQ_DDA)
%   incidence     incident field: |E|, phase and angular spectrum
%   apertures     extinction seen through central apertures of different sizes
%   emitted       "fraunhofer", "transmitted", "far_field", "quiver_ball"
%   coefficients  weight of the most excited eigenmodes vs wavelength (eig only)
%   coef_map      |c| of every labelled mode at chosen wavelengths: how strongly the
%                 incident field projects on each mode (c = V \ Einc); a (m_x, m_y)
%                 map, or |c| vs m for line arrays and cylinders (eig only)
%   optimal       best extinction over all unit-norm incident fields (any amplitude
%                 and phase on every dipole) vs the best single eigenmode, and their
%                 ratio. 1 means exciting one mode is optimal; > 1 means a mixture of
%                 non-orthogonal modes does better. Rebuilds A at n_points
%                 wavelengths (plus the coef_map ones), so it takes a while for
%                 large arrays; see theory.pdf, "Biorthogonal modes"
%   eig_band      extinction of the (mx,0) and (0,my) modes vs energy (eig only);
%                 only the direction(s) with more than one particle
% All wavelengths in the options are in nm. Mode labels are 0-based (mx,my);
% for line arrays and cylinders (m,0) and (0,m) both mean mode m along the line.

%% Load
if isstruct(file)
    S = file;
else
    path = char(file);
    if ~isfile(path), path = fullfile(dda_results_folder("real"), path); end
    S = load(path);
end
Input = S.Input;  C = S.C;  Array = S.Array;
opt = plotting_defaults();
if isfield(Input,'plotting'), opt = merge_deep(opt, Input.plotting); end
if nargin > 1 && ~isempty(plotting), opt = merge_deep(opt, plotting); end

D = prepare(Input, C, double(S.P), double(S.Einc), Array);
figs = struct();
negative_extinction_report(D);
dda_figure('docking', opt.docked);     % tabs in one window (Input.plotting.docked)

%% Plots
if opt.dipoles.enable
    figs.dipoles = dipole_overview(D, D.P, exp(-C.ext), 'driven response', opt.dipoles);
    if D.has_modes
        lins = [];
        for r = 1:size(opt.dipoles.modes,1)
            lin = mode_index(D, opt.dipoles.modes(r,:));
            if ~isnan(lin), lins(end+1) = lin; end %#ok<AGROW>
        end
        if opt.dipoles.n_strongest > 0
            lins = [lins, strongest_modes(D, opt.dipoles.n_strongest)];
        end
        for lin = unique(lins, 'stable')
            Pm = double(D.P_eigs{lin});   % stored single: |p|^2 ~ 1e-62 underflows
            Tm = exp(-D.extsall(:,lin).');
            figs.(sprintf('mode_%d',lin)) = dipole_overview(D, Pm, Tm, ...
                ['eigenmode ' mode_label(D,lin)], opt.dipoles);
        end
    elseif ~isempty(opt.dipoles.modes) && Input.inversion_tool == "eig"
        fprintf('No per-mode dipole moments in this file (C.P_eigs); mode plots skipped.\n');
    end
end

if opt.Q.enable
    wl0 = opt.Q.wl;
    if isempty(wl0), wl0 = D.lam(D.ipeak); end
    figs.Q = dda_figure([100 100 700 500], 'Name', 'Q fit');
    tlq = tiledlayout(figs.Q, 'flow');              % its own figure: never full
    tic
    try
        calc_Q_varQ_DDA(D.lam, C.ext, wl0, opt.Q.guess, tlq);
    catch err
        fprintf('Q fit failed (%s); continuing with the other plots.\n', err.message);
    end
    fprintf('Q fit took %.1f s\n', toc);
end

if opt.incidence.enable, figs.incidence = plot_incidence(D, opt.incidence.wl); end
if opt.apertures.enable, figs.apertures = plot_apertures(D, opt.apertures.sizes); end

if opt.emitted.enable
    wls = opt.emitted.wls;
    if isempty(wls), wls = D.lam(D.ipeak); end
    for wl = wls(:).'
        for t = string(opt.emitted.types)
            f = plot_emitted(D, t, wl, opt.emitted);
            if ~isempty(f), figs.(sprintf('emitted_%s_%dnm', t, round(wl))) = f; end
        end
    end
end

if opt.coefficients.enable
    if D.has_coefs, figs.coefficients = plot_coefficients(D, opt.coefficients.n_modes);
    else, fprintf('No eigenmode coefficients in this file; coefficient plot skipped.\n'); end
end
if opt.coef_map.enable
    if D.has_coefs, figs.coef_map = plot_coef_map(D, opt.coef_map.wls);
    else, fprintf('No eigenmode coefficients in this file; coefficient map skipped.\n'); end
end
if opt.optimal.enable
    wls = opt.optimal.wls;
    if isempty(wls) && opt.coef_map.enable, wls = opt.coef_map.wls; end
    figs.optimal = plot_optimal_excitation(D, opt.optimal.n_points, wls);
end
if opt.eig_band.enable
    if D.has_band, figs.eig_band = plot_eig_band(D);
    else, fprintf('No tracked mode extinctions in this file; band plot skipped.\n'); end
end

%% Export
if strlength(string(opt.export_folder)) > 0
    if ~isfolder(opt.export_folder), mkdir(opt.export_folder); end
    names = fieldnames(figs);
    for n = 1:numel(names)
        exportgraphics(figs.(names{n}), fullfile(opt.export_folder, [names{n} '.png']), ...
            'Resolution', 300, 'BackgroundColor', 'none');
    end
end
end

%% ========================================================================
%  Defaults (the same values DDA_Master sets; used for files without them)
%% ========================================================================
function o = plotting_defaults()
o.enable        = true;
o.docked        = true;
o.export_folder = "";
o.dipoles      = struct('enable',true, 'modes',[0 0; 1 0; 0 1; 0 2; 0 3; 0 4], ...
                        'n_strongest',0, 'profile_wls',1060.4, 'map_wl',[], 'surf_wls',[]);
o.dipoles.position_map = struct('quantity',"rss", 'column',0, 'spectral_unit',"nm", ...
    'spectral_lims',[], 'position_unit',"um", 'reverse_position',false, ...
    'spectral_axis',"auto", 'overlay',"none", 'overlay_color',[0.12 0.65 1], ...
    'colormap',"pink", 'clim',[], 'shading',"flat", 'colorbar_label',"", 'show_title',true);
o.Q            = struct('enable',true, 'wl',[], 'guess',30);
o.incidence    = struct('enable',true, 'wl',800);
o.apertures    = struct('enable',false, 'sizes',[1 10 20 30 40 60 80]);
o.emitted      = struct('enable',false, 'types',"fraunhofer", 'wls',630, ...
                        'angle_lim',2, 'fraunhofer_quantity',"P/alpha", ...
                        'anglex',linspace(-0.03,0.03,20), 'angley',[]);
o.coefficients = struct('enable',true, 'n_modes',10);
o.coef_map     = struct('enable',true, 'wls',[]);
o.optimal      = struct('enable',true, 'n_points',40, 'wls',[]);
o.eig_band     = struct('enable',true);
end

%% ========================================================================
%  Data preparation
%% ========================================================================
function D = prepare(Input, C, P, Einc, Array)
D.Input = Input;  D.C = C;  D.Array = Array;
D.Nk    = size(P,2);
D.N     = size(Array.r_i,1);
D.ncomp = size(P,1)/D.N;
D.P     = P;
D.Pp    = reshape(P,    D.ncomp, D.N, D.Nk);
D.Ep    = reshape(Einc, D.ncomp, D.N, D.Nk);
D.lam   = Input.range.wl(:).'*1e9;                % nm
D.eV    = 1239.84187./D.lam;
[~, D.ipeak] = max(abs(C.ext));

% Positions and grid indices (from the positions, so any ordering works)
r  = double(Array.r_i);
D.x = r(:,1)*1e6;  D.y = r(:,2)*1e6;              % um
D.z = zeros(D.N,1); if size(r,2) > 2, D.z = r(:,3)*1e6; end
tol = 1e-3*min(Input.p(1:2))*1e6;
[D.xv, ~, D.ix] = uniquetol_simple(D.x, tol);
[D.yv, ~, D.iy] = uniquetol_simple(D.y, tol);
D.is_grid = numel(D.xv)*numel(D.yv) == D.N && ...
            numel(unique(sub2ind([numel(D.yv) numel(D.xv)], D.iy, D.ix))) == D.N;

% Line arrays: the plots use only the direction along the line.
% Cylinders (along y) always form a line along x.
D.is_cyl = any(string(Input.dimensionality) == ["cyl_TM","cyl_TE"]);
if D.is_cyl || (numel(D.yv) == 1 && numel(D.xv) > 1)
    D.line = "x";  D.lv = D.xv;  D.il = D.ix;  D.lpos = D.x;
elseif numel(D.xv) == 1 && numel(D.yv) > 1
    D.line = "y";  D.lv = D.yv;  D.il = D.iy;  D.lpos = D.y;
else
    D.line = "";   D.lv = [];    D.il = [];    D.lpos = [];
end

% Mode data (eigen solver with tracking)
D.has_modes = isfield(C,'P_eigs');
D.has_band  = isfield(C,'extsall');
D.has_coefs = isfield(C,'coefs_ordered');
if isfield(C,'track')
    D.shape = C.track.shape;  D.lab_x = C.track.labels_x;  D.lab_y = C.track.labels_y;
elseif D.has_modes
    D.shape = size(C.P_eigs); D.lab_x = 0:D.shape(1)-1;     D.lab_y = 0:D.shape(2)-1;
elseif D.has_band
    sz = size(C.extsall);     D.shape = sz(2:end);
    D.lab_x = 0:D.shape(1)-1; D.lab_y = 0:D.shape(2)-1;
end
if D.has_band,  D.extsall = reshape(C.extsall, D.Nk, []); end
if D.has_modes, D.P_eigs  = C.P_eigs(:);                  end
if D.has_coefs, D.coefs   = reshape(C.coefs_ordered, [], D.Nk); end
D.label_kind = "mode";
if isfield(C,'track') && isfield(C.track,'basis'), D.label_kind = string(C.track.basis); end
end

function [u, iu, ic] = uniquetol_simple(v, tol)
[u, iu, ic] = unique(round(v/tol)*tol);
end

%% ========================================================================
%  Transmission above 1 (negative extinction)
%% ========================================================================
function negative_extinction_report(D)
% exp(-ext) > 1 means ext < 0: the array gives energy to the field. With
% passive particles that cannot happen, for any incident field, because
%   ext ~ sum Im(p.E*) = p' * (Im G - Im(1/alpha)) * p >= 0
% when every particle satisfies the optical theorem,
%   Im(1/alpha) <= -k^3/(6 pi eps0 n^2)      point dipoles (SI)
%   Im(1/alpha) <= -k^2/(4 eps0 n^2)          cyl_TM (per length)
%   Im(1/alpha) <= -k^2/(8 eps0 n^2)          cyl_TE (per length)
% This lists what in the run breaks it.
ext = double(D.C.ext(:).');
neg = ext < -1e-6*max(abs(ext));                 % (round-off of dark modes is not a cause)
if D.has_band
    eb  = double(D.extsall);
    neg = neg | any(eb < -1e-6*max(abs(eb(:)),[],'omitnan'), 2).';
end
if ~any(neg), return; end
fprintf(['\nTransmission exp(-ext) > 1: the extinction is negative at %d of %d wavelengths ' ...
         '(driven minimum %.3g).\n'], nnz(neg), D.Nk, min(ext));
I = D.Input;  found = false;

[excess, rad] = passivity_excess(D);
if ~isempty(excess)
    bad = any(excess > 1e-3*rad, 1);             % 1 x Nk (alpha is stored single)
    if any(bad)
        found = true;
        fprintf(['  The particles scatter more than they extinguish at %d of %d wavelengths: ' ...
                 'Im(1/alpha) is above the radiation-reaction limit by up to %.0f%%.\n'], ...
                nnz(bad), D.Nk, 100*max(excess(:)./rad(:)));
        if any(strcmpi(string(I.sample.type), ["Lorentzian","lorentzian","lorentzian Au","nonidentical"]))
            corr = "none";
            if isfield(I.particle,'lorentz_correction'), corr = string(I.particle.lorentz_correction); end
            fprintf(['  Cause: Lorentzian particles with lorentz_correction = "%s" have no radiation ' ...
                     'reaction.\n  Fix:   Input.particle.lorentz_correction = "radiative" ' ...
                     '(or "MLWA") and run again.\n'], corr);
        elseif isfield(I,'Hermitian') && I.Hermitian
            fprintf('  Cause: Input.Hermitian = true keeps only Re(alpha). Fix: Input.Hermitian = false.\n');
        elseif isfield(I.particle,'cyl_alpha_legacy') && I.particle.cyl_alpha_legacy && D.is_cyl
            fprintf('  Cause: Input.particle.cyl_alpha_legacy = true (old cylinder polarizability).\n');
        elseif string(I.sample.type) == "glass_Mie" && isfield(I,'particle_gain') && I.particle_gain < 0
            fprintf('  Cause: Input.particle_gain < 0 gives the particles gain (intended?).\n');
        else
            fprintf('  Check the polarizability model of Input.sample.type = "%s".\n', I.sample.type);
        end
    end
end
if isfield(I,'gain') && I.gain ~= 0
    found = true; fprintf('  Input.gain = %g: gain is intended to give negative extinction.\n', I.gain);
end
if isfield(I,'mirror_A') && any(I.mirror_A) && ...
        ~(isfield(D.C,'solver_version') && D.C.solver_version >= 4)
    found = true;
    fprintf(['  Input.mirror_A = [%s] in a result from before the mirror fix: particles on a mirror ' ...
             'plane were counted twice,\n  which makes the array non-passive. Run it again ' ...
             '(DDA_Master recalculates such files).\n'], num2str(I.mirror_A));
end
if isfield(I,'periodic_boundary') && any(I.periodic_boundary)
    found = true;
    if isfield(I,'periodic_cell') && ~isempty(I.periodic_cell)
        L = double(I.periodic_cell(:).').*[1 1];
        L = L(logical(I.periodic_boundary(:).'));
        if any(L > 0 & L < 1e-8)
            fprintf(['  Input.periodic_cell = [%s] m: far too small (nm typed as m?), so the ' ...
                     'periodic images sat on top of the particles.\n  The current ' ...
                     'Make_int_matrix_A_SI stops with an error for this.\n'], num2str(I.periodic_cell));
        end
    end
    fprintf(['  Input.periodic_boundary = [%s]: the truncated lattice sums are not energy ' ...
             'conserving; more image cells (Input.periodic_Ncells) reduce the error.\n'], ...
             num2str(I.periodic_boundary));
end
if D.is_cyl && isfield(I,'layers') && isfield(I.layers,'enable') && I.layers.enable
    nl = [I.layers.below_n(:); I.layers.above_n(:)];
    if any(imag(nl) < 0)
        found = true; fprintf('  A layer index has a negative imaginary part (gain).\n');
    end
end
if ~found
    fprintf('  No cause found in the inputs (polarizability, gain, mirrors, periodic boundary).\n');
end
fprintf('\n');
end

function [excess, rad] = passivity_excess(D)
% excess = Im(1/alpha) + radiation-reaction limit, components x types x Nk
% flattened to rows x Nk; > 0 means the particle is not passive.
excess = []; rad = [];
Ar = D.Array;
if isfield(Ar,'alpha_types'), A = Ar.alpha_types;
elseif isfield(Ar,'alpha'),   A = Ar.alpha;
else, return
end
if ~isfield(Ar,'k') || ~isfield(Ar,'n'), return; end
k = double(Ar.k(:).');  n = double(Ar.n(:).');
if any(imag(k) ~= 0), return; end                % complex frequencies: not a real-axis check
eps0 = 8.8541878188e-12;
dim  = string(D.Input.dimensionality);
switch dim
    case "cyl_TM", lim = k.^2./(4*eps0*n.^2);    comps = 1;
    case "cyl_TE", lim = k.^2./(8*eps0*n.^2);    comps = 1:2;
    case "1Dx",    lim = k.^3./(6*pi*eps0*n.^2); comps = 1;
    case "1Dy",    lim = k.^3./(6*pi*eps0*n.^2); comps = 2;
    case "1Dz",    lim = k.^3./(6*pi*eps0*n.^2); comps = 3;
    case "2D",     lim = k.^3./(6*pi*eps0*n.^2); comps = 1:2;
    otherwise,     lim = k.^3./(6*pi*eps0*n.^2); comps = 1:3;
end
if ~dda_alpha_is_SI(D.Input, Ar), lim = 2/3*k.^3; end   % Gaussian units
Nk = numel(k);
if ismatrix(A)                                   % cyl_TM: types x Nk
    inv_im = imag(1./double(A));
else
    inv_im = [];
    for c = comps(comps <= size(A,2))
        inv_im = [inv_im; imag(1./double(reshape(A(:,c,c,:), size(A,1), [])))]; %#ok<AGROW>
    end
end
if size(inv_im,2) ~= Nk, excess = []; return; end
rad    = repmat(lim, size(inv_im,1), 1);
excess = inv_im + rad;
end

%% ========================================================================
%  Dipole overview (old plot_dipole_moments and plot_dipoles)
%% ========================================================================
function fig = dipole_overview(D, Pm, T, what, o)
pabs = squeeze(sqrt(sum(abs(reshape(Pm, D.ncomp, D.N, D.Nk)).^2, 1)));   % N x Nk
if D.Nk == 1, pabs = pabs(:); end
[~, kmap] = min(abs(D.lam - pick(o.map_wl, D.lam(1))));

if D.line ~= ""                                 % line array or cylinders: two panels
    fig = dda_figure([100 100 1100 420], 'Name', what);
    tl  = tiledlayout(fig, 1, 2, 'TileSpacing','compact');
    nexttile(tl)
    transmission_plot(D, T)
    position_map(nexttile(tl), D, pabs, T, o.position_map, D.line)
    title(tl, sprintf('%s: %s', what, run_label(D)))
    surface_plots(D, Pm, pabs, what, o)
    return
end

fig = dda_figure([100 100 1200 600], 'Name', what);
tl  = tiledlayout(fig, 2, 3, 'TileSpacing','compact');

nexttile(tl)
transmission_plot(D, T)

position_map(nexttile(tl), D, pabs, T, o.position_map, "y")   % |p| per y row vs wavelength

nexttile(tl)                                     % summed over y, per x column
cols = sqrt(accum_rows(D.ix, pabs.^2));
imagesc(D.xv, D.lam, cols.'); set(gca,'YDir','normal')
xlabel('x (\mum)'); ylabel('Wavelength (nm)'); title('|p| summed over y'); colorbar; colormap(gca, pink)

nexttile(tl)                                     % profile along y through the centre column
[~, cx] = min(abs(D.xv - mean(D.xv)));
in = find(D.ix == cx);
[ys, o2] = sort(D.y(in)); in = in(o2);
wls = unique([D.lam(D.ipeak), o.profile_wls(:).'], 'stable');
hold on; leg = strings(0);
for wl = wls
    [~, kk] = min(abs(D.lam - wl));
    plot(ys, pabs(in, kk), 'LineWidth', 2);
    leg(end+1) = sprintf('%.1f nm', D.lam(kk)); %#ok<AGROW>
end
hold off; legend(leg); box on
xlabel('y (\mum)'); ylabel('|p|'); title('Profile through the centre column')

nexttile(tl)
map_plot(D, pabs(:, D.ipeak));
title(sprintf('|p| at the extinction peak, %.1f nm', D.lam(D.ipeak)))

nexttile(tl)
map_plot(D, pabs(:, kmap));
title(sprintf('|p| at %.1f nm', D.lam(kmap)))

title(tl, sprintf('%s: %s', what, run_label(D)))
surface_plots(D, Pm, pabs, what, o)
end

function position_map(ax, D, pabs, T, o, dir)
% |p| vs position along dir ("x" or "y") and wavelength or energy: the second
% panel of the dipole figures. Options (Input.plotting.dipoles.position_map):
%   quantity       "rss": sqrt(sum |p|^2) over the other direction, "mean": mean |p|,
%                  "center": the centre column (row), "column": column (row) no. o.column
%   column         0-based index across the line, for quantity "column"
%   spectral_unit  "nm" or "eV";  spectral_lims [] or [from to] in that unit
%   position_unit  "um" or "index" (particle number, 0-based)
%   reverse_position  true: position axis increasing downwards (left to right
%                  when it is the horizontal axis): 0 at the top
%   spectral_axis  "auto", "horizontal" or "vertical"
%   overlay        "none", "extinction" or "transmission" on a right-hand axis
%                  (needs the spectral axis horizontal; "auto" then gives that)
%   overlay_color  RGB;  colormap  name ("pink", "hot", ...) or an m x 3 matrix
%   clim           [] or [lo hi];  shading  "flat" or "interp" (smooth)
%   colorbar_label "" = automatic;  show_title  true/false
if dir == "y", idx = D.iy; oth = D.ix; pos = D.yv; ovals = D.xv; odir = "x";
else,          idx = D.ix; oth = D.iy; pos = D.xv; ovals = D.yv; odir = "y"; end
np = numel(pos);
q  = lower(string(o.quantity));
switch q
    case "rss"
        M = sqrt(rows_of(idx, pabs.^2, np));
        lab = '|p|';  if numel(ovals) > 1, lab = sprintf('|p| summed over %s', odir); end
    case "mean"
        M = rows_of(idx, pabs, np)./max(rows_of(idx, ones(D.N,1), np), 1);
        lab = '|p|';  if numel(ovals) > 1, lab = sprintf('mean |p| over %s', odir); end
    case {"center","centre","column"}
        if q == "column", c = min(max(round(o.column) + 1, 1), numel(ovals));
        else, [~, c] = min(abs(ovals - mean(ovals))); end
        sel = oth == c;
        M = sqrt(rows_of(idx(sel), pabs(sel,:).^2, np));
        lab = '|p|';
        if numel(ovals) > 1, lab = sprintf('|p| of %s column %d', odir, c - 1); end
    otherwise
        error('position_map.quantity must be "rss", "mean", "center" or "column", not "%s".', q);
end
M(rows_of(idx, ones(D.N,1), np) == 0, :) = NaN;     % positions without particles

if lower(string(o.position_unit)) == "index"
    pv = 0:np-1;  plab = sprintf('%s (NPs)', dir);
else
    pv = pos(:).';  plab = sprintf('%s (\\mum)', dir);
end
if lower(string(o.spectral_unit)) == "ev", sv = D.eV;  slab = 'Frequency (eV)';
else,                                      sv = D.lam; slab = 'Wavelength (nm)'; end
ov  = lower(string(o.overlay));
sa  = lower(string(o.spectral_axis));
hor = sa == "horizontal" || (sa == "auto" && (dir == "y" || ov ~= "none"));

if hor, draw_map(ax, sv, pv, M,   o.shading); xlabel(ax, slab); ylabel(ax, plab);
else,   draw_map(ax, pv, sv, M.', o.shading); xlabel(ax, plab); ylabel(ax, slab); end
if ~isempty(o.spectral_lims)
    if hor, xlim(ax, sort(o.spectral_lims)); else, ylim(ax, sort(o.spectral_lims)); end
end
if logical(o.reverse_position)
    if hor, set(ax, 'YDir', 'reverse'); else, set(ax, 'XDir', 'reverse'); end
end
cm = o.colormap;
if ~isnumeric(cm), cm = feval(char(cm), 256); end
colormap(ax, cm)
if ~isempty(o.clim), caxis(ax, o.clim); end
cb = colorbar(ax);
cbl = char(string(o.colorbar_label));  if isempty(cbl), cbl = lab; end
cb.Label.String = cbl;
if logical(o.show_title), title(ax, lab); end

if ov ~= "none"
    if ~hor
        fprintf('position_map: the %s overlay needs spectral_axis "horizontal"; skipped.\n', ov);
        return
    end
    if ov == "extinction", y = -log(T(:).'); ylb = 'Extinction';
    else,                  y = T(:).';       ylb = 'Transmission'; end
    oc = o.overlay_color;
    yyaxis(ax, 'right')
    plot(ax, sv, y, '-', 'Color', oc, 'LineWidth', 2)
    ylabel(ax, ylb);  ax.YAxis(2).Color = oc;
    yl = [min(0, min(y)) max(y)];  if diff(yl) <= 0, yl = yl + [-1 1]; end
    ylim(ax, yl + [0 0.05*diff(yl)])
    yyaxis(ax, 'left');  ax.YAxis(1).Color = [0.15 0.15 0.15];
end
end

function M = rows_of(idx, vals, n)
% sum the rows of vals that share idx, as an n x Nk array
M = zeros(n, size(vals,2));
if isempty(idx), return; end
for k = 1:size(vals,2), M(:,k) = accumarray(idx(:), vals(:,k), [n 1]); end
end

function draw_map(ax, x, y, C, shading)
% C(numel(y), numel(x)) at the points (x, y), which need not be evenly spaced
% (energy axis). "flat": one cell per point; "interp": smooth between points.
[x, ox] = sort(x(:).');  [y, oy] = sort(y(:).');  C = C(oy, ox);
if lower(string(shading)) == "interp" && numel(x) > 1 && numel(y) > 1
    surface(ax, x, y, zeros(size(C)), C, 'EdgeColor', 'none', 'FaceColor', 'interp');
    xl = x([1 end]);  yl = y([1 end]);
else
    xe = cell_edges(x);  ye = cell_edges(y);
    Cp = NaN(numel(ye), numel(xe));  Cp(1:end-1, 1:end-1) = C;
    surface(ax, xe, ye, zeros(size(Cp)), Cp, 'EdgeColor', 'none', 'FaceColor', 'flat');
    xl = xe([1 end]);  yl = ye([1 end]);
end
view(ax, 2);  box(ax, 'on');  set(ax, 'Layer', 'top', 'TickDir', 'out')
xlim(ax, xl);  ylim(ax, yl)
end

function e = cell_edges(c)
if numel(c) == 1, e = c + [-0.5 0.5]; return; end
m = (c(1:end-1) + c(2:end))/2;
e = [c(1) - (m(1) - c(1)), m, c(end) + (c(end) - m(end))];
end

function surface_plots(D, Pm, pabs, what, o)
if isempty(o.surf_wls), return; end              % old Input.plotdipoles
f2 = figure('Name', [what ' surfaces']);
t2 = tiledlayout(f2, 'flow');
Pr = reshape(Pm, D.ncomp, D.N, D.Nk);
for wl = o.surf_wls(:).'
    [~, kk] = min(abs(D.lam - wl));
    [~, c]  = max(sum(abs(Pr(:,:,kk)).^2, 2));
    nexttile(t2)
    surf_plot(D, pabs(:,kk), angle(Pr(c,:,kk)).');
    title(sprintf('%.1f nm', D.lam(kk)))
end
end

function transmission_plot(D, T)
plot(D.lam, T, 'LineWidth', 2)
xlabel('Wavelength (nm)'); ylabel('Transmission exp(-ext)'); box on
if any(T > 1)
    title('ext < 0 at some wavelengths: see the console', 'FontWeight','normal', 'Color',[0.8 0 0])
end
end

function s = run_label(D)
I = D.Input;
if D.is_cyl
    s = sprintf('%d cylinders, %s', D.N, I.dimensionality);
else
    s = sprintf('%d x %d, %s', I.arrayN(1), I.arrayN(min(2,end)), I.dimensionality);
end
if isfield(I,'gain') && I.gain ~= 0, s = sprintf('%s, gain %g', s, I.gain); end
end

function M = accum_rows(idx, vals)
% Sum the rows of vals (N x Nk) that share the same idx
M = zeros(max(idx), size(vals,2));
for k = 1:size(vals,2), M(:,k) = accumarray(idx, vals(:,k)); end
end

function map_plot(D, v)
if D.line ~= ""                                  % line: value vs position
    [ps, o] = sort(D.lpos);
    plot(ps, v(o), '.-', 'LineWidth', 1.5); box on; axis tight
    xlabel(sprintf('%s (\\mum)', D.line))
    return
end
if D.is_grid
    M = NaN(numel(D.yv), numel(D.xv));
    M(sub2ind(size(M), D.iy, D.ix)) = v;
    imagesc(D.xv, D.yv, M, 'Interpolation','bilinear'); set(gca,'YDir','normal')
else
    scatter(D.x, D.y, 36, v, 'filled'); axis equal tight
end
xlabel('x (\mum)'); ylabel('y (\mum)'); colorbar
end

function surf_plot(D, mag, ph)
if D.line ~= ""                                  % line: |p| vs position, coloured by phase
    scatter(D.lpos, mag, 36, ph, 'filled'); box on
    colormap(gca, hsv); clim([-pi pi]); colorbar
    xlabel(sprintf('%s (\\mum)', D.line)); ylabel('|p|')
    return
end
if D.is_grid
    M = NaN(numel(D.yv), numel(D.xv)); F = M;
    M(sub2ind(size(M), D.iy, D.ix)) = mag;
    F(sub2ind(size(F), D.iy, D.ix)) = ph;
    surf(D.xv, D.yv, M, F, 'EdgeColor','none')
else
    scatter3(D.x, D.y, mag, 36, ph, 'filled')
end
colormap(gca, hsv); clim([-pi pi]); colorbar
xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('|p|')
end

function v = pick(v, default)
if isempty(v), v = default; end
end

%% ========================================================================
%  Eigenmode helpers
%% ========================================================================
function lin = mode_index(D, row)
% row = [mx my] or [mx my pol] (0-based labels, pol 1..3); NaN if not available
lin = NaN;
if ~D.has_modes, return; end
if D.line ~= "" && any(row(1:2) == 0)            % line array: (m,0) and (0,m) = mode m along the line
    m = sum(row(1:2));
    if D.line == "x", row(1:2) = [m 0]; else, row(1:2) = [0 m]; end
end
ix = find(D.lab_x == row(1), 1);  iy = find(D.lab_y == row(2), 1);
if isempty(ix) || isempty(iy)
    fprintf('Mode (%d,%d) is not in this array; skipped.\n', row(1), row(2)); return
end
sub = ones(1, numel(D.shape)); sub(1) = ix; sub(2) = iy;
if numel(D.shape) >= 3
    if numel(row) >= 3
        sub(3) = row(3);
    elseif D.has_band                  % polarization with the strongest extinction
        best = -inf;
        for p = 1:D.shape(3)
            sub(3) = p;
            e = max(D.extsall(:, sub2ind_vec(D.shape, sub)));
            if e > best, best = e; bp = p; end
        end
        sub(3) = bp;
    end
end
lin = sub2ind_vec(D.shape, sub);
end

function lin = sub2ind_vec(shape, sub)
c = num2cell(sub); lin = sub2ind([shape 1], c{:});
end

function lins = strongest_modes(D, n)
if ~D.has_coefs, lins = []; return; end
a = abs(D.coefs(:, D.ipeak)); a(isnan(a)) = -inf;   % unassigned slots last
[~, lins] = maxk(a, n);
lins = lins(:).';
end

function s = mode_label(D, lin)
sub = cell(1, numel(D.shape)); [sub{:}] = ind2sub([D.shape 1], lin);
s = sprintf('(%d,%d)', D.lab_x(sub{1}), D.lab_y(sub{2}));
if numel(D.shape) >= 3 && D.shape(3) > 1, s = [s, ' ', comp_name(sub{3}), '-pol']; end
if numel(D.shape) >= 4 && D.shape(4) > 1, s = sprintf('%s, branch %d', s, sub{4}); end
end

%% ========================================================================
%  Incident field (old plot_incidence)
%% ========================================================================
function fig = plot_incidence(D, wl)
[~, kk] = min(abs(D.lam - wl));
E = D.Ep(:,:,kk);
[~, c] = max(sum(abs(E).^2, 2));
v = E(c,:).';
fig = dda_figure([100 100 1200 400], 'Name','incident field');
tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact');

nexttile(tl); map_plot(D, abs(v));
title(sprintf('|E_%s|', comp_name(c))); clim([0 max(abs(v))])
nexttile(tl); map_plot(D, angle(v)/pi); colormap(gca, hsv); clim([-1 1]); title('Phase/\pi')

nexttile(tl)
if D.line ~= ""                                  % line: spectrum along the line only
    [th, Pw] = angular_spectrum_line(D, v, kk);
    plot(th, Pw, 'LineWidth', 1.5); box on; axis tight
    xlabel(sprintf('\\theta_%s (deg)', D.line)); ylabel('|FFT|^2'); title('Angular spectrum')
elseif D.is_grid
    [TX, TY, Pw] = angular_spectrum(D, v, kk);
    surf(TX, TY, Pw, 'EdgeColor','none'); view(2); axis tight
    xlabel('\theta_x (deg)'); ylabel('\theta_y (deg)'); title('Angular spectrum'); colorbar
else
    axis off; text(0.1, 0.5, 'Angular spectrum needs a rectangular grid')
end
title(tl, sprintf('%.1f nm, \\theta = %.1f, \\phi = %.1f, \\psi = %.1f', ...
    D.lam(kk), D.Input.theta, D.Input.phi, D.Input.psi))
end

function [th, Pw] = angular_spectrum_line(D, v, kk, nn)
% |FFT|^2 along a line array, on propagating angles only
if nargin < 4, nn = 1024; end
L = zeros(numel(D.lv), 1);
L(D.il) = v;                                     % (stacked particles: the last one)
d = mean(diff(D.lv))*1e-6;                       % spacing along the line [m]
if isempty(d) || ~isfinite(d) || d <= 0, d = D.Input.p(1 + (D.line == "y")); end
F  = fftshift(fft(L, nn));
k0 = 2*pi/(D.lam(kk)*1e-9);
kl = 2*pi*(-nn/2:nn/2-1).'/(nn*d);
ok = abs(kl) < k0;
th = asind(kl(ok)/k0);
Pw = abs(F(ok)).^2/nn;
end

function [TX, TY, Pw, F] = angular_spectrum(D, v, kk, nn)
% |FFT|^2 of a grid map, on propagating angles only
if nargin < 4, nn = 1024; end
M = zeros(numel(D.yv), numel(D.xv));
M(sub2ind(size(M), D.iy, D.ix)) = v;
F  = fftshift(fft2(M, nn, nn));
k0 = 2*pi/(D.lam(kk)*1e-9);
fx = (-nn/2:nn/2-1)/(nn*D.Input.p(1));  kx = 2*pi*fx;
fy = (-nn/2:nn/2-1)/(nn*D.Input.p(2));  ky = 2*pi*fy;
okx = abs(kx) < k0;  oky = abs(ky) < k0;
[TX, TY] = meshgrid(asind(kx(okx)/k0), asind(ky(oky)/k0));
F  = F(oky, okx);
Pw = abs(F).^2/nn;
end

%% ========================================================================
%  Apertures (old plot_transmission_apertures)
%% ========================================================================
function fig = plot_apertures(D, sizes)
cx = (max(D.ix)+1)/2;  cy = (max(D.iy)+1)/2;
fig = figure('Name','apertures'); hold on
leg = strings(0);
for a = sizes(:).'
    sel = abs(D.ix - cx) < a/2 & abs(D.iy - cy) < a/2;
    if ~any(sel) || all(sel), continue; end
    plot(D.lam, dda_extinction(D.P, reshape(D.Ep,[],D.Nk), D.Input, D.Array, sel), 'LineWidth', 1.5)
    leg(end+1) = sprintf('%d x %d particles (%.1f x %.1f \\mum)', min(a,max(D.ix)), ...
        min(a,max(D.iy)), min(a,max(D.ix))*D.Input.p(1)*1e6, min(a,max(D.iy))*D.Input.p(2)*1e6); %#ok<AGROW>
end
plot(D.lam, D.C.ext, 'k', 'LineWidth', 2); leg(end+1) = "whole array";
hold off; legend(leg, 'Location','best'); box on
xlabel('Wavelength (nm)'); ylabel('Extinction')
title('Extinction through central apertures')
end

%% ========================================================================
%  Emitted / transmitted fields (old plot_emitted)
%% ========================================================================
function fig = plot_emitted(D, type, wl, o)
fig = [];
[~, kk] = min(abs(D.lam - wl));
[~, c]  = max(sum(abs(D.Ep(:,:,kk)).^2, 2));         % dominant incident component
switch type
    case {"fraunhofer","transmitted"}
        if ~D.is_grid
            fprintf('"%s" needs a rectangular grid; skipped.\n', type); return
        end
        p   = D.Pp(c,:,kk).';
        alp = local_alpha(D, c, kk);
        if type == "fraunhofer"
            if string(o.fraunhofer_quantity) == "P", q = p; else, q = p./alp; end
            [TX, TY, ~, F] = angular_spectrum(D, q, kk);
            fig = figure('Name', sprintf('Fraunhofer %.0f nm', D.lam(kk)));
            tl = tiledlayout(fig, 1, 2);
            nexttile(tl); surf(TX,TY,abs(F),'EdgeColor','none'); view(2); colorbar; title('|F|')
            lim_axes(o.angle_lim)
            nexttile(tl); surf(TX,TY,real(F),'EdgeColor','none'); view(2); colorbar; title('Re F')
            lim_axes(o.angle_lim)
            title(tl, sprintf('Fraunhofer pattern of %s, %.1f nm', o.fraunhofer_quantity, D.lam(kk)))
        else
            [TX, TY, ~, Finc] = angular_spectrum(D, D.Ep(c,:,kk).', kk);
            [~,  ~,  ~, Femi] = angular_spectrum(D, p./alp/4/pi, kk);
            Ft = Finc - Femi;
            fig = figure('Name', sprintf('transmitted %.0f nm', D.lam(kk)));
            tl = tiledlayout(fig, 1, 2);
            nexttile(tl); surf(TX,TY,abs(Ft),'EdgeColor','none'); view(2); colorbar
            title('|E_{transmitted}|'); lim_axes(o.angle_lim)
            nexttile(tl); surf(TX,TY,abs(Ft).^2./abs(Finc).^2,'EdgeColor','none'); view(2)
            colorbar; clim([0 1]); title('Transmission'); lim_axes(o.angle_lim)
            title(tl, sprintf('Transmitted field, %.1f nm', D.lam(kk)))
        end

    case "far_field"
        ax = o.anglex(:).';  ay = o.angley;  if isempty(ay), ay = ax.'; end
        [th, ph] = anglexy_to_theta_phi(ax, ay(:));
        rh = reshape(thetaphi_to_rhat(th, ph), [], 3);
        E  = far_field(D, rh, kk);
        fig = figure('Name', sprintf('far field %.0f nm', D.lam(kk)));
        imagesc(rad2deg(ax), rad2deg(ay), reshape(sqrt(sum(abs(E).^2,2)), numel(ay), numel(ax)), ...
            'Interpolation','bilinear'); set(gca,'YDir','normal'); colorbar
        xlabel('\theta_x (deg)'); ylabel('\theta_y (deg)')
        title(sprintf('|E| emitted, %.1f nm', D.lam(kk)))

    case "quiver_ball"
        [TH, PH] = meshgrid(linspace(0,pi,60), linspace(0,2*pi,60));
        rh = reshape(thetaphi_to_rhat(TH, PH), [], 3);
        E  = far_field(D, rh, kk);
        mag = reshape(sqrt(sum(abs(E).^2,2)), size(TH));
        X = reshape(rh(:,1),size(TH)); Y = reshape(rh(:,2),size(TH)); Z = reshape(rh(:,3),size(TH));
        fig = figure('Name', sprintf('quiver ball %.0f nm', D.lam(kk)));
        surf(X, Y, Z, mag/max(mag(:)), 'EdgeColor','none', 'FaceColor','interp'); hold on
        s = 1:6:numel(TH);  Er = real(E(s,:)); Er = Er./max(vecnorm(Er,2,2), eps);
        quiver3(X(s).', Y(s).', Z(s).', Er(:,1).', Er(:,2).', Er(:,3).', 0.3, 'k')
        hold off; axis equal tight; view(-45,20); colorbar
        xlabel('x'); ylabel('y'); zlabel('z')
        title(sprintf('Emitted far field, %.1f nm, %d particles', D.lam(kk), D.N))

    otherwise
        error('Unknown emitted plot type "%s".', type);
end
end

function E = far_field(D, rh, kk)
% E(rhat) ~ sum_j exp(-i k rhat.r_j) (rhat x p_j) x rhat = Q - rhat (rhat.Q)
p3 = dipoles_3D(D, kk);                                   % N x 3
r  = [D.x D.y D.z]*1e-6;
k  = double(D.Array.k(kk));
Q  = exp(-1i*k*(rh*r.'))*p3;                              % M x 3
E  = Q - rh.*sum(rh.*Q, 2);
end

function p3 = dipoles_3D(D, kk)
p  = D.Pp(:,:,kk).';                                      % N x ncomp
p3 = zeros(D.N, 3);
switch string(D.Input.dimensionality)
    case "3D",     p3 = p;
    case "2D",     p3(:,1:2) = p;
    case "1Dx",    p3(:,1) = p;
    case "1Dy",    p3(:,2) = p;
    case "1Dz",    p3(:,3) = p;
    case "cyl_TM", p3(:,2) = p;        % assumed: field along the cylinder axis (y)
    case "cyl_TE", p3(:,[1 3]) = p;    % line dipoles in the xz plane
end
end

function alp = local_alpha(D, c, kk)
% Polarizability of each particle for component c (Array.alpha is N x 3 x 3 x Nk)
comp = c;
switch string(D.Input.dimensionality)
    case "1Dx", comp = 1;  case "1Dy", comp = 2;  case "1Dz", comp = 3;
end
alp = double(squeeze(dda_alpha(D.Array, (1:D.N).', kk, comp, comp)));   % full or compact storage
end

function lim_axes(lim)
xlim([-lim lim]); ylim([-lim lim]); xlabel('\theta_x (deg)'); ylabel('\theta_y (deg)')
end

%% ========================================================================
%  Eigenmode coefficients and band diagram (old plot_coefs, eig_banddiagram)
%% ========================================================================
function fig = plot_coef_map(D, wls)
% |c| of the labelled modes, c = V \ Einc: how much of the incident field each
% mode takes. One tile per wavelength (default: the extinction peak). The
% polarizations and branches of a label are combined as sqrt(sum |c|^2).
% Each tile gives the strongest mode and its share |c_max|^2 / sum |c|^2
% (100 % = no other mode is excited). The eigenvectors are not orthogonal, so
% the shares are a measure of selectivity, not exact energy fractions.
if isempty(wls), wls = D.lam(D.ipeak); end
A  = reshape(abs(D.coefs).^2, [D.shape(1:2), prod(D.shape(3:end)), D.Nk]);
A  = sqrt(reshape(sum(A, 3, 'omitnan'), D.shape(1), D.shape(2), D.Nk));   % nx x ny x Nk
[lx, ox] = sort(D.lab_x);  [ly, oy] = sort(D.lab_y);   % Fourier labels may be signed
A  = A(ox, oy, :);
show = [D.shape(1) > 1, D.shape(2) > 1];
nt = numel(wls);
fig = dda_figure([100 100 min(1500, 470*nt) 380], 'Name', 'mode coefficients');
tl  = tiledlayout(fig, 1, nt, 'TileSpacing', 'compact');
for wl = wls(:).'
    [~, kk] = min(abs(D.lam - wl));
    ax = nexttile(tl);
    if all(show)                                   % 2D: map over (m_x, m_y)
        imagesc(ax, lx, ly, A(:,:,kk).');  set(ax, 'YDir', 'normal');  axis(ax, 'tight')
        xlabel(ax, 'Mode number m_x');  ylabel(ax, 'Mode number m_y')
        colormap(ax, flipud(pink));  c = colorbar(ax);  c.Label.String = '|c|';
    else                                           % line array: |c| vs m
        if show(2), m = ly; a = squeeze(A(1,:,kk)); d = 'y'; else, m = lx; a = squeeze(A(:,1,kk)); d = 'x'; end
        bar(ax, m, a, 0.8, 'FaceColor', [0.6 0.35 0.3], 'EdgeColor', 'none');  box(ax, 'on')
        xlabel(ax, sprintf('Mode number m_%s', d));  ylabel(ax, '|c|');  axis(ax, 'tight')
    end
    Ak = A(:,:,kk);  P = Ak(:).^2;  P(isnan(P)) = 0;
    [pmax, l] = max(P);  [ix, iy] = ind2sub(size(Ak), l);
    if all(show) || ~any(show), lab = sprintf('(%d,%d)', lx(ix), ly(iy));
    elseif show(2),            lab = sprintf('%d', ly(iy));
    else,                      lab = sprintf('%d', lx(ix));
    end
    title(ax, sprintf('%.1f nm: mode %s takes %.1f %%', D.lam(kk), lab, 100*pmax/max(sum(P), eps)))
end
title(tl, sprintf('Incident field %s on the modes (%s labels): %s', field_label(D.Input), ...
      D.label_kind, run_label(D)))
end

function s = field_label(I)
% polarization, angles and profile of the incident field, e.g.
% "linear (0, 0, 90)°, Gaussian w0 = 80 µm, HG01"
s = sprintf('%s (%g, %g, %g)°', char(string(I.polarization)), I.theta, I.phi, I.psi);
if isfield(I, 'preset_field') && any(strcmpi(char(string(I.preset_field)), {'Gaussian', 'gaussian_proper'}))
    s = sprintf('%s, Gaussian w0 = %.3g µm', s, I.Gaussian_w0*1e6);
    if I.Gaussian_z ~= 0, s = sprintf('%s, z = %.3g µm', s, I.Gaussian_z*1e6); end
    s = sprintf('%s, HG%d%d', s, I.Hermite_order_x, I.Hermite_order_y);
else
    s = sprintf('%s, plane wave', s);
end
end

%% ========================================================================
%  Optimal excitation vs single-mode excitation
%% ========================================================================
function fig = plot_optimal_excitation(D, n_points, wls)
% For a unit-norm incident field E (one value per dipole component) the
% extinction is  ext = K * Im(E' * A^-1 * E),  K = 4*pi*k/cell  (/(4*pi*eps0*n^2)
% for SI polarizabilities): the normalisation of dda_extinction, in which
% area*I0 = cell*|E|^2. Hence
%   one eigenmode v_a (|v_a| = 1):  ext_a   = K * Im(1/lambda_a)   (= C.extsall)
%   best of all fields:             ext_opt = K * max eig of H,
%                                   H = (A^-1 - A^-H)/(2i)  (Hermitian)
% ext_opt >= max_a ext_a. Equality holds when the eigenvectors are orthogonal
% (A normal); A is complex symmetric, not normal, so a mixture of modes can
% do better. The overlap |u' v_best| between the optimal field u and the best
% mode says how close the two are (1 = the same field). At the peak of an
% isolated mode (Re 1/lambda = 0) the ratio is (1 + sqrt(K))/2, K = Petermann
% factor |v|^4/|v.'v|^2 of the mode; drawn dotted for comparison.
I = D.Input;  Ar = D.Array;  Nk = D.Nk;
Kf = 4*pi*double(Ar.k(:).');
if D.is_cyl, cell_area = I.p(1); else, cell_area = I.p(1)*I.p(2); end
Kf = Kf/cell_area;
if dda_alpha_is_SI(I, Ar), Kf = Kf/(4*pi*8.8541878188e-12)./double(Ar.n(:).').^2; end

has_vals = isfield(D.C,'values') && ~isempty(D.C.values) && ~all(isnan(D.C.values(:)));
single_all = NaN(1, Nk);
if has_vals
    single_all = Kf.*max(imag(1./double(D.C.values)), [], 1);
end

% Wavelengths at which A is rebuilt
ks = unique(round(linspace(1, Nk, min(n_points, Nk))));
[~, kp] = max(single_all);  if has_vals, ks = [ks kp]; end
ks = [ks D.ipeak];
for wl = wls(:).', [~, kk] = min(abs(D.lam - wl)); ks = [ks kk]; end %#ok<AGROW>
ks = unique(ks);
nks = numel(ks);
opt_ext = NaN(1,nks);  sgl_ext = NaN(1,nks);  overlap = NaN(1,nks);  Kp = NaN(1,nks);
fprintf('Optimal excitation: %d wavelengths, %d unknowns ... ', nks, D.N*D.ncomp);  t0 = tic;
ws = warning('off', 'MATLAB:nearlySingularMatrix');  restore = onCleanup(@() warning(ws));
wsel = ks == D.ipeak;                              % wavelengths reported in the console
if has_vals, wsel = wsel | ks == kp; end
for wl = wls(:).', [~, kk] = min(abs(D.lam - wl)); wsel = wsel | ks == kk; end
for j = 1:nks
    kk = ks(j);
    A  = double(build_A(I, Ar, kk));
    n  = size(A,1);
    Ai = A \ eye(n);
    H  = (Ai - Ai')/2i;  H = (H + H')/2;
    [U, h] = eig(H, 'vector');
    [hmax, i] = max(real(h));
    u = U(:,i);
    opt_ext(j) = Kf(kk)*hmax;
    if has_vals                                    % best mode by inverse iteration
        lam = double(D.C.values(:,kk));
        [~, a] = max(imag(1./lam));  la = lam(a);
        v = u;
        for it = 1:3
            v = (A - la*(1 + 1e-12)*eye(n)) \ v;  v = v/norm(v);
        end
    else                                           % no saved eigenvalues: solve here
        [V, lam] = eig(A, 'vector');
        [~, a] = max(imag(1./lam));  la = lam(a);  v = V(:,a)/norm(V(:,a));
    end
    sgl_ext(j) = Kf(kk)*imag(1/la);
    overlap(j) = abs(u'*v);
    Kp(j)      = 1/abs(v.'*v)^2;                   % Petermann factor of the best mode
end
fprintf('%.1f s\n', toc(t0));
if ~has_vals, single_all(ks) = sgl_ext; end
ratio = opt_ext./sgl_ext;

fig = dda_figure([100 100 1000 420], 'Name', 'optimal excitation');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl);  hold(ax, 'on');  box(ax, 'on')
plot(ax, D.lam, D.C.ext, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'DisplayName', 'this incident field')
if has_vals
    plot(ax, D.lam, single_all, 'LineWidth', 1.5, 'DisplayName', 'best single mode')
else
    plot(ax, D.lam(ks), sgl_ext, '-', 'LineWidth', 1.5, 'DisplayName', 'best single mode')
end
plot(ax, D.lam(ks), opt_ext, 'o-', 'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', 'optimal field')
hold(ax, 'off');  axis(ax, 'tight');  legend(ax, 'Location', 'best')
xlabel(ax, 'Wavelength (nm)');  ylabel(ax, 'Extinction (unit-norm field)')
ax = nexttile(tl);
yyaxis(ax, 'left');  plot(ax, D.lam(ks), ratio, 'o-', 'MarkerSize', 3, 'LineWidth', 1.5);  hold(ax, 'on')
plot(ax, D.lam(ks), (1 + sqrt(Kp))/2, ':', 'LineWidth', 1.2);  hold(ax, 'off')
ylabel(ax, 'ext_{opt} / ext_{best mode}');  yl = ylim(ax);  ylim(ax, [min(1, yl(1)) max(yl(2), 1.05)])
yyaxis(ax, 'right'); plot(ax, D.lam(ks), overlap, '.--', 'LineWidth', 1)
ylabel(ax, '|u^H v_{best}|');  ylim(ax, [0 1.02])
legend(ax, {'optimal / best mode', '(1+\surdK)/2 of that mode', 'overlap'}, 'Location', 'best')
xlabel(ax, 'Wavelength (nm)');  box(ax, 'on');  xlim(ax, [min(D.lam) max(D.lam)])
[~, jp] = max(opt_ext);
title(tl, sprintf(['Optimal vs single-mode excitation: at the optimum peak (%.1f nm) ' ...
      'a mode mixture gives %.3g x the best mode, overlap %.2f. %s'], ...
      D.lam(ks(jp)), ratio(jp), overlap(jp), run_label(D)))
for j = find(wsel)
    fprintf('  %8.2f nm  %7.4f eV  ratio %.4f  overlap %.3f  Petermann K %.3f\n', ...
        D.lam(ks(j)), D.eV(ks(j)), ratio(j), overlap(j), Kp(j));
end
end

function A = build_A(Input, Array, kk)
% The solver's interaction matrix (as in Calculate_DDA_eigenmodes)
switch string(Input.dimensionality)
    case "cyl_TM", A = Make_cyl_TM_matrix_A_SI(Input, Array, kk);
    case "cyl_TE", A = Make_cyl_TE_matrix_A_SI(Input, Array, kk);
    otherwise,     A = Make_int_matrix_A_SI(Input, Array, kk);
end
end

function fig = plot_coefficients(D, n)
w = abs(D.coefs).^2 ./ sum(abs(D.coefs).^2, 1, 'omitnan');
lins = strongest_modes(D, n);
fig = figure('Name','eigenmode coefficients'); hold on
leg = strings(0);
for lin = lins
    plot(D.lam, w(lin,:), 'LineWidth', 1.5); leg(end+1) = mode_label(D, lin); %#ok<AGROW>
end
hold off; legend(leg, 'Location','best'); box on
xlabel('Wavelength (nm)'); ylabel('|c|^2 / \Sigma|c|^2')
title(sprintf('Most excited modes at %.1f nm', D.lam(D.ipeak)))
end

function fig = plot_eig_band(D)
% Extinction of the (mx,0) and (0,my) modes vs energy; only the directions
% with more than one particle (one panel for line arrays and cylinders)
E3 = sum(reshape(D.extsall, D.Nk, D.shape(1), D.shape(2), []), 4, 'omitnan');   % sum pol/branch
show = [D.shape(1) > 1, D.shape(2) > 1];
if ~any(show), show(1) = true; end               % a single particle
fig = dda_figure([100 100 550*nnz(show) 450], 'Name','eigenmode bands');
tl = tiledlayout(fig, 1, nnz(show), 'TileSpacing','compact');
if show(1)
    nexttile(tl)
    imagesc(D.lab_x, D.eV, E3(:,:,1)); set(gca,'YDir','normal')
    xlabel('Mode number m_x'); ylabel('Energy (eV)'); c = colorbar; c.Label.String = 'Extinction';
    if D.is_cyl, title('Modes m_x'); else, title('Modes (m_x, 0)'); end
end
if show(2)
    nexttile(tl)
    imagesc(D.lab_y, D.eV, reshape(E3(:,1,:), D.Nk, [])); set(gca,'YDir','normal')
    xlabel('Mode number m_y'); ylabel('Energy (eV)'); title('Modes (0, m_y)'); c = colorbar; c.Label.String = 'Extinction';
end
title(tl, sprintf('Eigenmode bands: %s', run_label(D)))
end

%% ========================================================================
function s = comp_name(c)
names = 'xyz'; s = names(min(c,3));
end

function a = merge_deep(a, b)
f = fieldnames(b);
for n = 1:numel(f)
    if isfield(a, f{n}) && isstruct(a.(f{n})) && isstruct(b.(f{n}))
        a.(f{n}) = merge_deep(a.(f{n}), b.(f{n}));
    else
        a.(f{n}) = b.(f{n});
    end
end
end
