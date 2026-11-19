function out = dda_mode_basis(varargin)
% DDA_MODE_BASIS  Mode labels (mx, my, polarization) from basis functions at
% the particle positions. Shared by Track_eigenmodes (real frequencies) and
% complex_freqs_plots_from_file, so both give a mode the same label.
%
%   M = dda_mode_basis(r, nx, ny, basis, Input)    the basis
%       r      N x 2 or N x 3 particle positions
%       nx,ny  number of labels in x and y (Input.arrayN(1:2))
%       basis  "auto" (Fourier along periodic directions, sine otherwise) |
%              "sine" | "fourier" | "k" (plane-wave coupling, below)
%       M.B (N x nx*ny, unit columns, mx fastest), M.lab_x, M.lab_y, M.used,
%       M.nx, M.ny, M.kind ("basis" or "k")
%
%   S = dda_mode_basis(M, V, ncomp)                basis scores
%       V      Ndof x K vectors, components interleaved per particle
%       S      nx*ny*ncomp x K: fraction of each vector's power on each label
%              (rows: component blocks of nx*ny labels; columns sum to 1)
%
% Basis functions: sin((m+1)*pi*u) along open directions (the standing-wave /
% DST-I basis, m = number of nodes) and exp(2i*pi*f*x/L) along periodic
% directions. They are evaluated at the actual positions, so shifted,
% half-shifted or non-square lattices work too.
%
% "k" labels (plane-wave coupling): the coupling of a vector v to the plane
% waves exp(i(+-kx x +- ky y)), C(|kx|,|ky|) = sum over the signs of
% |v.' e_k|^2, is computed on a |k| grid. Each label m owns the k values
% nearest to k_m, the coupling peak of the standing wave sin((m+1) pi u) on
% the same positions (calibration), and its score is the largest coupling in
% that interval; the label with the largest score is therefore the k label of
% complex_freqs_plots_from_file. Input.track.k_step = true first undoes the pi
% phase steps of half-shifted lattices (the step-matched, HG-shaped field).
% Along periodic directions with more than one cell "k" falls back to "auto".
%
%   K = dda_mode_basis("k_setup", r, [nx ny], os)   k grids and calibration
%       os: grid points per pi/(L+2d) (default 40 for line arrays, 10 otherwise)
%   C = dda_mode_basis("k_spectrum", K, V)          nky x nkx x size(V,2) coupling

if isstruct(varargin{1})
    out = basis_score(varargin{:});
elseif (ischar(varargin{1}) || isstring(varargin{1})) && ...
        any(strcmp(char(varargin{1}), {'k_setup', 'k_spectrum'}))
    if strcmp(char(varargin{1}), 'k_setup'), out = k_setup(varargin{2:end});
    else,                                    out = k_spectrum(varargin{2:end});
    end
else
    out = position_basis(varargin{:});
end
end

%% ------------------------------------------------------------------------
function M = position_basis(r, nx, ny, basis, Input)
r = double(r);
pbc = [false false];
if isfield(Input,'periodic_boundary'), pbc = logical(Input.periodic_boundary(:).'); end
basis = lower(string(basis));
if basis == "k"
    if any(pbc & [nx ny] > 1)
        fprintf('dda_mode_basis: "k" labels are not defined along periodic directions; using "auto".\n');
        basis = "auto";
    else
        M = k_basis(r, nx, ny, Input);
        return
    end
end
switch basis
    case "auto",    per = pbc;
    case "sine",    per = [false false];
    case "fourier", per = [true true];
    otherwise, error('dda_mode_basis:basis','The basis must be "auto", "sine", "fourier" or "k".');
end
[Bx, lab_x] = dir_basis(r(:,1), nx, per(1));
[By, lab_y] = dir_basis(r(:,2), ny, per(2));
B = reshape(Bx .* permute(By,[1 3 2]), size(r,1), nx*ny);   % mx fastest, as sub2ind([nx ny])
nrm = vecnorm(B); nrm(nrm == 0) = 1;
names = ["sine","fourier"];
M = struct('B', B ./ nrm, 'lab_x', lab_x, 'lab_y', lab_y, ...
           'used', names(per(1)+1) + "/" + names(per(2)+1), 'nx', nx, 'ny', ny, 'kind', "basis");
end

function [b, labels] = dir_basis(x, n, periodic)
labels = 0:n-1;
if n <= 1
    b = ones(numel(x),1);
    return
end
xmin = min(x); span = max(x) - xmin;
if span <= 0
    b = [ones(numel(x),1), zeros(numel(x), n-1)];
    return
end
d = span/(n-1);                                % lattice spacing from the positions
if periodic
    labels = [0:ceil(n/2)-1, -floor(n/2):-1];  % spatial frequency in units of 1/L
    b = exp(2i*pi*(x - xmin).*labels/(n*d));
else
    u = (x - xmin + d)/(span + 2*d);           % DST-I positions for a regular grid
    b = sin(pi*u.*(1:n));                      % column m+1: m nodes
end
end

%% ------------------------------------------------------------------------
function score = basis_score(M, V, ncomp)
if isfield(M, 'kind') && strcmp(char(M.kind), 'k')
    score = k_score(M, V, ncomp);
    return
end
B = M.B;
ncell = size(B,2);
score = zeros(ncell*ncomp, size(V,2));
for c = 1:ncomp
    Vc = V(c:ncomp:end, :);                    % components are interleaved per particle
    score((c-1)*ncell + (1:ncell), :) = abs(B'*Vc).^2;
end
tot = sum(score,1);
tot(tot == 0) = 1;
score = score ./ tot;
end

%% ------------------------------------------------------------------------
%  Plane-wave coupling ("k") labels
function M = k_basis(r, nx, ny, Input)
ks = false;
if isfield(Input, 'track') && isstruct(Input.track) && isfield(Input.track, 'k_step')
    ks = logical(Input.track.k_step);
end
% k grid: as fine as the complex plots for line arrays (identical labels); for
% 2D arrays 4 points per pi/(L+2d) (the plots use 10), because every
% eigenvector of every frequency is scored
os = 4;  if nx == 1 || ny == 1, os = 40; end
M = struct('kind', "k", 'K', k_setup(r, [nx ny], os), 'k_step', ks, ...
           'lab_x', 0:nx-1, 'lab_y', 0:ny-1, 'nx', nx, 'ny', ny, 'B', []);
M.used = "k";
if ks, M.used = "k (steps undone)"; end
end

function K = k_setup(r, n, os)
r = double(r);
if nargin < 3 || isempty(os), os = 10; if any(n == 1), os = 40; end, end
N = size(r, 1);
for dd = 1:2
    c = r(:,dd);  span = max(c) - min(c);
    if n(dd) <= 1 || span <= 0
        K.k{dd} = 0;  K.km{dd} = 0;  K.dk(dd) = Inf;  K.E{dd} = ones(N,1);  K.grp{dd} = 1;
        continue
    end
    d  = span/(n(dd)-1);
    kg = (0:pi/(span+2*d)/os:1.05*pi/d).';
    Eg = exp(1i*c*kg.');                           % N x nk
    u  = (c - min(c) + d)/(span + 2*d);
    Sk = abs(Eg.'*sin(pi*u*(1:n(dd)))).^2;         % real standing waves: +k and -k equal
    [~, im] = max(Sk, [], 1);
    K.k{dd} = kg;  K.km{dd} = kg(im);  K.E{dd} = Eg;
    K.dk(dd) = 0.8*pi/(span + 2*d);                % half-width of a coupling peak
    [~, K.grp{dd}] = min(abs(kg - K.km{dd}(:).'), [], 2);   % label owning each k
end
% grid of the particles (separable, fast spectrum)
tol = 1e-9;
[xu, ~, K.ix] = unique(round(r(:,1)/tol)*tol);
[yu, ~, K.iy] = unique(round(r(:,2)/tol)*tol);
K.isgrid = numel(xu)*numel(yu) == N && numel(unique(K.iy + (K.ix-1)*numel(yu))) == N;
K.size = [numel(yu) numel(xu)];
K.Exu = exp(1i*xu(:)*K.k{1}(:).');               % nxu x nkx
K.Eyu = exp(1i*yu(:)*K.k{2}(:).');               % nyu x nky
% sign flips of the half-period shifts (array_shifts): a gap of more than
% 1.25 median spacings starts a new segment
K.flip = ones(N, 1);
for dd = 1:2
    if n(dd) <= 1, continue; end
    cu = unique(round(r(:,dd)/tol)*tol);
    if numel(cu) < 3, continue; end
    gaps  = diff(cu);
    stepc = cu([false; gaps > 1.25*median(gaps)]);
    seg   = sum(r(:,dd) >= stepc(:).' - tol, 2);
    K.flip = K.flip .* (-1).^seg;
end
end

function C = k_spectrum(K, V)
% coupling |v.' e_k|^2 summed over the four sign combinations, nky x nkx x nv
nv = size(V, 2);
C  = zeros(numel(K.k{2}), numel(K.k{1}), nv);
if K.isgrid
    Ex = K.Exu;  Ey = K.Eyu;
    G  = zeros(K.size);
    idx = sub2ind(K.size, K.iy, K.ix);
    for j = 1:nv
        G(idx) = V(:,j);
        A1 = Ey.'*G;  A2 = Ey'*G;
        C(:,:,j) = abs(A1*Ex).^2 + abs(A1*conj(Ex)).^2 + abs(A2*Ex).^2 + abs(A2*conj(Ex)).^2;
    end
else
    Ex = K.E{1};  Ey = K.E{2};
    for j = 1:nv
        v = V(:,j);
        C(:,:,j) = abs(Ey.'*(v.*Ex)).^2 + abs(Ey.'*(v.*conj(Ex))).^2 + ...
                   abs(Ey'*(v.*Ex)).^2  + abs(Ey'*(v.*conj(Ex))).^2;
    end
end
end

function score = k_score(M, V, ncomp)
K = M.K;  nx = M.nx;  ny = M.ny;  ncell = nx*ny;
score = zeros(ncell*ncomp, size(V,2));
for c = 1:ncomp
    Vc = V(c:ncomp:end, :);
    if M.k_step, Vc = Vc .* K.flip; end
    for j0 = 1:200:size(V,2)                      % in blocks: the spectra can be large
        cols = j0:min(size(V,2), j0+199);
        C  = k_spectrum(K, Vc(:,cols));           % nky x nkx x nc
        nc = numel(cols);
        Cx = zeros(size(C,1), nx, nc);
        for m = 1:nx
            sel = K.grp{1} == m;
            if any(sel), Cx(:,m,:) = max(C(:,sel,:), [], 2); end
        end
        S = zeros(ny, nx, nc);
        for m = 1:ny
            sel = K.grp{2} == m;
            if any(sel), S(m,:,:) = max(Cx(sel,:,:), [], 1); end
        end
        score((c-1)*ncell + (1:ncell), cols) = reshape(permute(S, [2 1 3]), ncell, nc);
    end
end
tot = sum(score, 1);
tot(tot == 0) = 1;
score = score ./ tot;
end
