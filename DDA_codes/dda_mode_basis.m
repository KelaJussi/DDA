function out = dda_mode_basis(varargin)
% DDA_MODE_BASIS  Mode labels (mx, my, polarization) from basis functions at
% the particle positions. Shared by Track_eigenmodes (real frequencies) and
% complex_freqs_plots_from_file, so both give a mode the same label.
%
%   M = dda_mode_basis(r, nx, ny, basis, Input)    the basis
%       r      N x 2 or N x 3 particle positions
%       nx,ny  number of labels in x and y (Input.arrayN(1:2))
%       basis  "auto" (Fourier along periodic directions, sine otherwise) |
%              "sine" | "fourier"
%       M.B (N x nx*ny, unit columns, mx fastest), M.lab_x, M.lab_y, M.used,
%       M.nx, M.ny
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

if isstruct(varargin{1})
    out = basis_score(varargin{:});
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
switch basis
    case "auto",    per = pbc;
    case "sine",    per = [false false];
    case "fourier", per = [true true];
    otherwise, error('dda_mode_basis:basis','The basis must be "auto", "sine" or "fourier".');
end
[Bx, lab_x] = dir_basis(r(:,1), nx, per(1));
[By, lab_y] = dir_basis(r(:,2), ny, per(2));
B = reshape(Bx .* permute(By,[1 3 2]), size(r,1), nx*ny);   % mx fastest, as sub2ind([nx ny])
nrm = vecnorm(B); nrm(nrm == 0) = 1;
names = ["sine","fourier"];
M = struct('B', B ./ nrm, 'lab_x', lab_x, 'lab_y', lab_y, ...
           'used', names(per(1)+1) + "/" + names(per(2)+1), 'nx', nx, 'ny', ny);
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
