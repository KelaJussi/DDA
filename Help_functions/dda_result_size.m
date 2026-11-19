function R = dda_result_size(Input)
% DDA_RESULT_SIZE  Estimated size of an eigen-solver run, in bytes.
% Used by Calculate_DDA_eigenmodes (its own check) and dda_size_check (before
% submitting), so both use the same numbers.
%
%   R.ndof, R.nk       matrix size and number of frequencies
%   R.lean             complex-frequency run that keeps only what the plots use
%   R.track            mode tracking on
%   R.results          what the run returns and saves (C, P, incident field)
%   R.vectors          eigenvectors held in memory while tracking (0 otherwise)
%   R.eig_workspace    memory one worker needs for one frequency (eig of an Ndof x Ndof matrix)
switch string(Input.dimensionality)
    case "3D", nc = 3;
    case {"2D","cyl_TE"}, nc = 2;
    otherwise, nc = 1;
end
ndof = nc*Input.N;
nk   = numel(Input.range.real_parts)*numel(Input.range.imag_parts);
if isfield(Input.range, 'k_N'), nk = Input.range.k_N; end
is_complex = any(Input.range.imag_parts ~= 0);

e = "auto";
if isfield(Input, 'track') && isfield(Input.track, 'enable'), e = Input.track.enable; end
if isstring(e) || ischar(e), track = ~is_complex; else, track = logical(e); end
complex_all = isfield(Input, 'save') && isfield(Input.save, 'complex_all') && Input.save.complex_all;
lean = is_complex && ~track && ~complex_all;

if lean
    results = 8*ndof*nk + 40*nk;                  % min eigenvector (single), min eigenvalue, degeneracy
else
    results = 64*ndof*nk;                         % all eigenvalues, min eigenvectors, coefficients,
end                                               % mode extinctions, P and incident field
vectors = 0;
if track
    results = results + 56*ndof*nk;               % ordered eigenvalues, coefficients, extinctions, labels
    store = ~isfield(Input, 'track') || ~isfield(Input.track, 'store_P_eigs') || Input.track.store_P_eigs;
    if store, results = results + 8*ndof^2*nk; end   % per-mode dipole moments (single)
    vectors = 16*ndof^2*nk;                       % all eigenvectors, while tracking
end
R = struct('ndof', ndof, 'nk', nk, 'lean', lean, 'track', track, ...
           'results', results, 'vectors', vectors, 'eig_workspace', 6*16*ndof^2);
end
