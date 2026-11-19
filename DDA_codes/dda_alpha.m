function a = dda_alpha(Array, ii, kk, ci, cj)
% DDA_ALPHA  Polarizability of particles ii at frequency index kk.
%
%   a = dda_alpha(Array, ii, kk)           N x 3 x 3 tensors (or N x 1 for cyl_TM)
%   a = dda_alpha(Array, ii, kk, ci, cj)   only components ci, cj
%
% Reads both storage forms:
%   full     Array.alpha        N x 3 x 3 x Nk   (or N x Nk for cyl_TM)
%   compact  Array.alpha_types  M x 3 x 3 x Nk   one tensor per particle type,
%            Array.alpha_type   N x 1            the type of each particle
% Identical particles share one type, so the compact form is N times smaller
% (dda_compact_alpha builds it, dda_expand_alpha restores the full form).
if isfield(Array, 'alpha_types')
    A = Array.alpha_types;  t = Array.alpha_type(ii);
else
    A = Array.alpha;        t = ii;
end
if ismatrix(A)                               % cyl_TM: scalar polarizability
    a = A(t, kk);
elseif nargin < 4
    a = A(t, :, :, kk);
else
    a = A(t, ci, cj, kk);
end
end
