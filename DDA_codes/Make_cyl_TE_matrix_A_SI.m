function [A]=Make_cyl_TE_matrix_A_SI(Input,Array,kk)
% MAKE_CYL_TE_MATRIX_A_SI  Interaction matrix of infinite cylinders (along y),
% field in the plane perpendicular to them (cyl_TE), 2N x 2N.
%
% Each cylinder is a line dipole with the two components (p_x, p_z); the
% unknowns are interleaved per cylinder, [p_x1 p_z1 p_x2 p_z2 ...], x along the
% array and z along the layer normal. With E = k^2/(eps0 n^2) * G * p:
%
%   A(i,i) = inv(alpha_i) - k^2/(eps0 n^2) * Gr(r_i, r_i)
%   A(i,j) =              - k^2/(eps0 n^2) * (G0(r_i - r_j) + Gr(r_i, r_j))
%
%   G0 = (I + grad grad/k^2) i/4*H0(kR)
%      = i/4*[ (H0 - H1/(kR)) I + H2 * Rhat Rhat ]        (2 x 2, in the xz plane)
%
% Gr is what a layered background adds (Input.layers, see dda_layers and
% dda_layer_green); it is zero without layers.

kh   = double(Array.k(kk));
nh   = double(Array.n(kk));
eps0 = 8.8541878188e-12; %F/m
N    = Input.N;

r_i = double(Array.r_i);
if any(abs(r_i(:,2) - r_i(1,2)) > 1e-12)
    error('Make_cyl_TE_matrix_A_SI:y', ['cyl_TE: the cylinders lie along y, so all particles must ' ...
        'have the same y. Let the array run along x (arrayN = [Nx 1 1]).']);
end
dx = r_i(:,1) - r_i(:,1).';
dz = r_i(:,3) - r_i(:,3).';
R  = sqrt(dx.^2 + dz.^2);
off = R > 1e-14;
if nnz(~off) > N
    error('Make_cyl_TE_matrix_A_SI:overlap', 'Two cylinders are at the same position.');
end
R(~off) = 1/abs(kh);                           % placeholder, removed below
u  = kh*R;
H0 = besselh(0,1,u);  H1 = besselh(1,1,u);  H2 = besselh(2,1,u);
cx = dx./R;  cz = dz./R;
iso = 1i/4*(H0 - H1./u);
Gxx = (iso + 1i/4*H2.*cx.*cx).*off;
Gzz = (iso + 1i/4*H2.*cz.*cz).*off;
Gxz = (1i/4*H2.*cx.*cz).*off;
Gzx = Gxz;

% Layered background: reflected field, including the one a cylinder gets back itself
L = dda_layers(Input, Array);
if L.enable
    Gr  = dda_layer_green(kh/nh, nh, L, dx, repmat(L.z, 1, N), repmat(L.z.', N, 1), 'TE');
    Gxx = Gxx + Gr.xx;  Gzz = Gzz + Gr.zz;  Gxz = Gxz + Gr.xz;  Gzx = Gzx + Gr.zx;
end

c = -kh^2/(eps0*nh^2);
A = zeros(2*N, 2*N);
A(1:2:end, 1:2:end) = c*Gxx;
A(1:2:end, 2:2:end) = c*Gxz;
A(2:2:end, 1:2:end) = c*Gzx;
A(2:2:end, 2:2:end) = c*Gzz;

% Self term: inverse of the 2 x 2 polarizability of each cylinder
a = double(dda_alpha(Array, (1:N).', kk));     % N x 2 x 2
a = reshape(a, N, 2, 2);
for n = 1:N
    ii = 2*n + [-1 0];
    A(ii, ii) = A(ii, ii) + inv(reshape(a(n,:,:), 2, 2));
end
end
