function F = dda_layer_field(k0, nh, L, kx, z, pol)
% DDA_LAYER_FIELD  Field of a plane wave inside the host layer of a layered background.
%
%   F = dda_layer_field(k0, nh, L, kx, z, pol)
%
% A plane wave of unit electric-field amplitude comes down from the top
% half-space of the stack L (see dda_layers) with in-plane wavenumber kx.
% Returns the field at the heights z (above the lower interface of the host
% layer), without the common factor exp(i*kx*x):
%   pol = 'TM'  F.y        E_y / E_inc                      (s polarization)
%   pol = 'TE'  F.x, F.z   E_x, E_z / E_inc                 (p polarization)
% For TE the sign is such that without interfaces E = -(kz, kx)/k * exp(-i*kz*z):
% the incident field is E_inc * (-cos(theta), -sin(theta)) for kx = k sin(theta).
%
% Transfer matrices for u (E_y or H_y) and w = (kz/eta)(up - down), eta = 1 (TM)
% or n^2 (TE), both continuous across the interfaces; the solution is started
% in the bottom half-space with a purely outgoing wave.
is_te = strcmpi(pol, 'TE');
z  = double(z(:));
kzf = @(n) pick(sqrt((n*k0)^2 - kx^2));
eta = @(n) 1 + is_te*(n^2 - 1);

nb = L.nb;  na = L.na;
kb = kzf(nb(end));
u = 1;  w = -kb/eta(nb(end));                        % outgoing (downward) wave in the substrate
for j = numel(nb)-1:-1:1                             % up through the layers below
    [u, w] = step(u, w, kzf(nb(j)), eta(nb(j)), L.db(j));
end
kh = kzf(nh);  eh = eta(nh);
uz = zeros(size(z));  wz = uz;
for m = 1:numel(z)                                   % inside the host layer
    [uz(m), wz(m)] = step(u, w, kh, eh, z(m));
end
[u, w] = step(u, w, kh, eh, L.d);
for j = 1:numel(na)-1                                % up through the layers above
    [u, w] = step(u, w, kzf(na(j)), eta(na(j)), L.da(j));
end
nt = na(end);  kt = kzf(nt);
a_inc = (u - w*eta(nt)/kt)/2;                        % incoming (downward) amplitude at the top
if is_te
    F.x = wz*nt/(k0*a_inc);                          % E_x = w/(omega eps0), E_inc = H_inc k_t/(omega eps0 n_t^2)
    F.z = -kx*uz*nt/(nh^2*k0*a_inc);
else
    F.y = uz/a_inc;
end
end

function [u2, w2] = step(u, w, kz, eta, d)
c = cos(kz*d);  s = sin(kz*d);
u2 = u*c + 1i*eta/kz*w*s;
w2 = w*c + 1i*kz/eta*u*s;
end

function kz = pick(kz)
% outgoing or decaying branch
if abs(real(kz)) >= abs(imag(kz))
    if real(kz) < 0, kz = -kz; end
elseif imag(kz) < 0
    kz = -kz;
end
end
