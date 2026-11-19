function Array = Lorentzian_polarizability_DDA(Array, Input)
% LORENTZIAN_POLARIZABILITY_DDA  MJH Lorentzian scatterers in the solver's SI units.
% Used when Input.sample.type = "Lorentzian" (replaces Make_Lorentzian_scatterer
% and Make_nonidentical_scatterers).
%
% Model, per axis c = x, y, z (MJH):
%   alpha_c(w) = A0_c / (w0_c - w - i*gamma_c),  w0 = 2*pi*c0/lambda_peak,
%                                                gamma = 1/(2*pi*tau)
% with w = 2*pi*c/wl the angular frequency of Input.range.eV (complex
% frequencies work).
%
% Units. MJH's A0 values give a Gaussian (CGS-style) volume polarizability
% alpha_G [m^3] for the old solver, A = 1/alpha_G - k^2*G with
% G = exp(ikR)/R*[...] and no host factor. The SI solver uses
% A = 1/alpha_SI - k^2/(4*pi*eps0*n^2)*G, so the same physics needs
%   alpha_SI = 4*pi*eps0*n_h^2 * alpha_G      (p_SI = 4*pi*eps0*n_h^2 * p_G)
% n_h is the host index at each frequency (BK7 or Input.n, as in the Mie path).
% This is also the factor between the Gaussian radiation reaction
% (2i/3)k^3 and the SI one i*k^3/(6*pi*eps0*n^2) used in Make_int_matrix_A_SI.
%
% Inputs (Input.particle.*)
%   scatter_strength  A0: scalar, [Ax Ay Az] or N x 3 (one row per particle).
%                     Legacy forms: 2 values = A0 ramps linearly over the
%                     wavelengths (as before); 6 values = [~ Ax ~ Ay ~ Az] with
%                     MLWA (as before; prefer lorentz_correction = "MLWA").
%   tau               lifetime(s) [s]: scalar, [x y z] or N x 3
%   lambda_peak_vec   resonance wavelength(s) [m]: scalar, [x y z] or N x 3
%   lorentz_units     "gaussian" (default): A0 as in MJH's code, converted
%                     "SI": A0 already gives alpha_SI [C m^2 V^-1], no conversion
%   lorentz_correction "none" (default, MJH's plain Lorentzian)
%                     "radiative": 1/alpha -> 1/alpha - i k^3/(6 pi eps0 n^2),
%                                  energy-conserving (optical theorem)
%                     "MLWA":      also the dynamic depolarization -k^2/(4 pi eps0 n^2 a)
%   lorentz_mlwa_a    particle size a for MLWA [m] (default 100e-9, as before)
%
% Outputs: Array.alpha_types (types x 3 x 3 x Nk, SI, zero off-diagonal) and
% Array.alpha_type (N x 1), see dda_alpha; Array.n, Array.k = n*w/c,
% Array.alpha_units = "SI".

eps0 = 8.8541878188e-12;
c0   = 2.99792458e8;
if isfield(Array, 'radius')
    warning('Lorentzian_polarizability_DDA:radius', 'Edited radii have no effect on Lorentzian particles.');
end
N    = Input.N;
P    = Input.particle;
units      = lower(string(getopt(P, 'lorentz_units', "gaussian")));
correction = lower(string(getopt(P, 'lorentz_correction', "none")));

if any(Input.dimensionality == ["cyl_TM","cyl_TE"])
    error('Lorentzian_polarizability_DDA:cyl', ...
        'Lorentzian particles are point dipoles; cyl_TM and cyl_TE need a 2D (per length) polarizability.');
end

%% Frequencies and host medium (same as Mie_polarizability_DDA)
eV = Input.range.eV(:).';
if isfield(Array,'wl')
    w = 2*pi*c0./Array.wl(:).';                      % as in MJH's code (wl = 1239.84187e-9./eV),
else                                                 % so his fitted A0/peaks/tau keep their meaning
    w = 2*pi*c0*eV/1239.84187e-9;
end
Nk = numel(w);
if Input.substrate == "BK7"
    n = sqrt(epsilon_BK7(eV));
else
    n = Input.n*ones(1, Nk);
end
Array.n = n;
Array.k = n.*w/c0;                                   % medium wavenumber n*2*pi/wl

%% Parameters per particle and axis
A0 = P.scatter_strength;
A0_ramp = ones(1, Nk);
if numel(A0) == 2                                     % legacy: linear ramp over wavelength
    A0_ramp = fliplr(linspace(A0(1), A0(2), Nk)) / A0(1);
    A0 = A0(1);
elseif numel(A0) == 6                                 % legacy: [~ Ax ~ Ay ~ Az] meant MLWA
    warning('Lorentzian_polarizability_DDA:legacy6', ...
        ['6-element scatter_strength is read as [~ Ax ~ Ay ~ Az] with MLWA, as in the old code. ' ...
         'Prefer 3 values and Input.particle.lorentz_correction = "MLWA".']);
    A0 = A0([2 4 6]);
    correction = "mlwa";
end
A0  = per_particle(A0,                N, 'scatter_strength');
tau = per_particle(P.tau,             N, 'tau');
lp  = per_particle(P.lambda_peak_vec, N, 'lambda_peak_vec');
% One polarizability per distinct parameter set (particle type)
[~, first, typ] = unique([A0 tau lp], 'rows');
A0 = A0(first,:);  tau = tau(first,:);  lp = lp(first,:);
Nt = numel(first);
w0    = 2*pi*c0./lp;                                  % Nt x 3
gamma = 1./(2*pi*tau);

%% Polarizability, Nt x 3 x Nk
alpha = complex(zeros(Nt, 3, Nk));
for c = 1:3
    alpha(:,c,:) = reshape(A0(:,c).*A0_ramp ./ (w0(:,c) - w - 1i*gamma(:,c)), Nt, 1, Nk);
end

nr = reshape(n,       1, 1, Nk);
kr = reshape(Array.k, 1, 1, Nk);
switch units
    case "gaussian", alpha = 4*pi*eps0*nr.^2.*alpha;
    case "si"        % already SI
    otherwise, error('Input.particle.lorentz_units must be "gaussian" or "SI".');
end

switch correction
    case "none"
    case "radiative"
        alpha = 1./(1./alpha - 1i*kr.^3./(6*pi*eps0*nr.^2));
    case "mlwa"
        a = getopt(P, 'lorentz_mlwa_a', 100e-9);
        alpha = 1./(1./alpha - 1i*kr.^3./(6*pi*eps0*nr.^2) - kr.^2./(4*pi*eps0*nr.^2*a));
    otherwise
        error('Input.particle.lorentz_correction must be "none", "radiative" or "MLWA".');
end

%% Passivity check at real frequencies: Im(1/alpha) <= -k^3/(6 pi eps0 n^2)
if correction == "none" && all(abs(imag(eV)) < eps)
    excess = imag(1./alpha) + real(kr).^3./(6*pi*eps0*nr.^2);
    frac = mean(excess(:) > 0);
    if frac > 0
        % printed, not warning(): DDA_Master switches warnings off
        fprintf(['Lorentzian_polarizability_DDA: at %.0f%% of the frequencies the Lorentzian scatters ' ...
                 'more than it extinguishes (no radiation reaction), so the extinction can be ' ...
                 'negative (transmission > 1).\n  Use Input.particle.lorentz_correction = "radiative".\n'], 100*frac);
    end
end

%% Store per type (compact form, read with dda_alpha): diagonal tensors,
%  exact zeros off the diagonal; alpha_type(i) is the type of particle i
Array.alpha_types = zeros(Nt, 3, 3, Nk, 'like', alpha);
for c = 1:3
    Array.alpha_types(:,c,c,:) = reshape(alpha(:,c,:), Nt, 1, 1, Nk);
end
Array.alpha_type = typ(:);
if isfield(Array, 'alpha'), Array = rmfield(Array, 'alpha'); end
Array.alpha_units     = "SI";
Array.lambda_peak_vec = P.lambda_peak_vec;
end

%% ------------------------------------------------------------------------
function V = per_particle(v, N, name)
% scalar -> all particles and axes; [x y z] -> all particles; N x 3 as is
if isscalar(v)
    V = repmat(double(v), N, 3);
elseif numel(v) == 3
    V = repmat(double(v(:).'), N, 1);
elseif isequal(size(v), [N 3])
    V = double(v);
else
    error('Lorentzian_polarizability_DDA:size', ...
        'Input.particle.%s must be a scalar, 1 x 3 or %d x 3 (one row per particle).', name, N);
end
end

function v = getopt(s, name, default)
if isfield(s, name), v = s.(name); else, v = default; end
end
