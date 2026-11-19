function [A]=Make_cyl_TM_matrix_A_SI(Input,Array,kk)
% MAKE_CYL_TM_MATRIX_A_SI  Interaction matrix of infinite cylinders (along y),
% field along the cylinders (cyl_TM), N x N.
%
%   A(i,i) = 1/alpha_i - k^2/(eps0 n^2) * Gr(r_i, r_i)
%   A(i,j) =           - k^2/(eps0 n^2) * (G0(r_i - r_j) + Gr(r_i, r_j))
%
% G0 = i/4*H0(kR) is the 2D Green's function of the host medium. Gr is what a
% layered background adds (Input.layers, see dda_layers and dda_layer_green);
% it is zero without layers.

kh   = double(Array.k(kk));
nh   = double(Array.n(kk));
eps0 = 8.8541878188e-12; %F/m
N    = Input.N;

r_i = double(Array.r_i);
R = sqrt((r_i(:,1)-r_i(:,1).').^2 + (r_i(:,2)-r_i(:,2).').^2 + (r_i(:,3)-r_i(:,3).').^2);
R(R<1e-14)=0;

% 2D Green function of the host medium, no self-interaction
G = 1i/4 * besselh(0,1,kh*R);
G(1:N+1:end) = 0;

% Layered background: reflected field, including the one a cylinder gets back itself
L = dda_layers(Input, Array);
if L.enable
    dx = r_i(:,1) - r_i(:,1).';
    Gr = dda_layer_green(kh/nh, nh, L, dx, repmat(L.z, 1, N), repmat(L.z.', N, 1), 'TM');
    G  = G + Gr.yy;
end

A = -kh.^2 * G ./ (eps0.*nh.^2);

% Add self term (polarizability)
a_self = dda_alpha(Array, (1:N).', kk);   % scalar alpha (full or compact storage)
A(1:N+1:end) = A(1:N+1:end) + 1./reshape(double(a_self), 1, []);
end
