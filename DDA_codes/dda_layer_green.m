function G = dda_layer_green(k0, nh, L, dx, zi, zj, pol)
% DDA_LAYER_GREEN  Reflected part of the 2D Green's function in a layered medium.
%
%   G = dda_layer_green(k0, nh, L, dx, zi, zj, pol)
%
% Line sources along y inside the host layer of the stack L (see dda_layers).
% Returns only what the interfaces add to the homogeneous Green's function of
% the host medium, for source j and observation point i in the host layer:
%   dx = x_i - x_j,  zi, zj = heights above the lower interface of the host layer
% (arrays of the same size, one entry per pair; dx = 0 with zi = zj is the
% field a source gets back from the interfaces, the self term).
%
%   pol = 'TM'  field along the cylinders (cyl_TM):  G.yy
%   pol = 'TE'  field in the xz plane (cyl_TE):      G.xx, G.xz, G.zx, G.zz
%
% Same normalisation as the homogeneous part used by the matrix builders:
%   E = k^2/(eps0*n^2) * G * p,   k = nh*k0,
%   G0_yy = i/4*H0(kR),  G0 = (I + grad grad/k^2) i/4*H0(kR) in the plane.
%
% k0 is the vacuum wavenumber and may be complex (complex frequencies).
%
% Method. Plane-wave (Sommerfeld) integral over kx with the generalized
% reflection coefficients of the two stacks, s-polarized for TM and
% p-polarized for TE:
%   G = i/(4 pi) int dkx exp(i kx dx)/kz * [ Rd e^{ikz(zi+zj)} + Ru e^{ikz(2d-zi-zj)}
%                                          + Rd Ru e^{2ikz d} 2cos(kz(zi-zj)) ] / (1 - Rd Ru e^{2ikz d})
% (TE: with the derivative factors of the in-plane dyadic). The path runs
% below the real kx axis, deeper than every branch point and guided-mode pole,
% so it stays valid for complex frequencies (Im(w) < 0 moves the poles of
% lossless guided modes below the axis). Beyond all singularities the path
% turns into the complex plane, where the integrand decays exponentially.
% kz is followed continuously along the path from kz = k at kx = 0.

if ~L.enable, error('dda_layer_green:off', 'The layered background is not enabled.'); end
is_te = strcmpi(pol, 'TE');
sz  = size(dx);
dx  = double(dx(:));  zi = double(zi(:));  zj = double(zj(:));
k0  = double(k0);     nh = double(nh);
kh  = nh*k0;
d   = L.d;
nall = [nh, L.nb, L.na];
kall = nall*k0;

%% Integration path
dref  = max([max(abs(dx)); 5/abs(real(k0))]);
delta = 3/dref;                                   % clearance from the singularities
h0    = 1.25*max([0, -imag(kall)]) + delta;       % depth below the real axis
kre   = abs(real(kall));
t_r   = 0.5*max(min(kre), 0.1*abs(real(k0)));     % ramp from kx = 0 down to the depth h0
T1    = 1.5*max(abs(kall)) + 2*delta;             % beyond every branch point and guided pole
zs    = min([zi + zj; 2*d - zi - zj]);            % shortest way via an interface
S_end = max(60/zs, 4*T1);

[xg, wg]   = gauss_legendre(8);
[xg2, wg2] = gauss_legendre(12);

% core: kx = t - i*h(t), 0 <= t <= T1
edges = unique([linspace(0, t_r, max(2, ceil(t_r/delta)+1)), ...
                linspace(t_r, T1, max(2, ceil((T1-t_r)/delta)+1))]);
[t, wt] = panels(edges, xg, wg);
hh   = h0*min(t/t_r, 1);
kx_c = t - 1i*hh;
w_c  = wt.*(1 - 1i*h0/t_r*(t < t_r));

% tails: from kx = T1 - i*h0 along exp(+i*pi/4) and exp(-i*pi/4), graded panels
ne    = ceil(log(S_end/delta)/log(1.25)) + 1;
edges = [0, delta*1.25.^(0:ne)];
[s, ws] = panels(edges, xg2, wg2);
ph   = exp(1i*pi/4);
kx_p = T1 - 1i*h0 + s*ph;        w_p = ws*ph;
kx_m = T1 - 1i*h0 + s*conj(ph);  w_m = ws*conj(ph);

nc = numel(kx_c);  np = numel(kx_p);
kx = [kx_c; kx_p; kx_m];
w  = [w_c;  w_p;  w_m];

%% kz in every layer, continuous along the path
kz = zeros(numel(kx), numel(kall));
for m = 1:numel(kall)
    a = follow_sqrt(kall(m)^2 - kx_c.^2, kall(m));
    b = follow_sqrt(kall(m)^2 - kx_p.^2, a(end));
    c = follow_sqrt(kall(m)^2 - kx_m.^2, a(end));
    kz(:,m) = [a; b; c];
end
kzh = kz(:,1);
kzb = kz(:, 1 + (1:numel(L.nb)));
kza = kz(:, 1 + numel(L.nb) + (1:numel(L.na)));

%% Generalized reflection coefficients seen from the host layer
Rd = stack_reflection(kzh, nh, kzb, L.nb, L.db, is_te);
Ru = stack_reflection(kzh, nh, kza, L.na, L.da, is_te);
e2d = exp(2i*kzh*d);
M   = 1./(1 - Rd.*Ru.*e2d);

%% Phase factors exp(+-i kx |dx|) for the distinct |dx|
tol = 1e-13;
[adx, ~, idx] = unique(round(abs(dx)/tol)*tol);
adx = adx(:).';
Pe = [2*cos(kx_c*adx); exp(1i*kx_p*adx);  exp(-1i*kx_m*adx)];      % even in kx
if is_te
    Po = [2i*sin(kx_c*adx); exp(1i*kx_p*adx); -exp(-1i*kx_m*adx)]; % odd in kx
end

%% One spectral weight per distinct pair of heights
[zz, ~, iz] = unique(round([zi zj]/tol)*tol, 'rows');
pref = 1i/(4*pi)*w./kzh;
if is_te
    G = struct('xx', zeros(sz), 'xz', zeros(sz), 'zx', zeros(sz), 'zz', zeros(sz));
else
    G = struct('yy', zeros(sz));
end
sgn = sign(dx);
for m = 1:size(zz, 1)
    za = zz(m,1);  zb = zz(m,2);
    T1_ = M.*Rd.*exp(1i*kzh*(za + zb));
    T2_ = M.*Ru.*exp(1i*kzh*(2*d - za - zb));
    T3_ = M.*Rd.*Ru.*exp(1i*kzh*(2*d + za - zb));
    T4_ = M.*Rd.*Ru.*exp(1i*kzh*(2*d - za + zb));
    sel = find(iz == m);
    col = idx(sel);
    if is_te
        fxx = pref.*kzh.^2/kh^2.*(-T1_ - T2_ + T3_ + T4_);
        fzz = pref.*kx.^2/kh^2.*( T1_ + T2_ + T3_ + T4_);
        fxz = -pref.*kx.*kzh/kh^2.*(T1_ - T2_ + T3_ - T4_);
        fzx =  pref.*kx.*kzh/kh^2.*(T1_ - T2_ - T3_ + T4_);
        gxx = fxx.'*Pe;  gzz = fzz.'*Pe;  gxz = fxz.'*Po;  gzx = fzx.'*Po;
        G.xx(sel) = gxx(col);
        G.zz(sel) = gzz(col);
        G.xz(sel) = gxz(col).'.*sgn(sel);
        G.zx(sel) = gzx(col).'.*sgn(sel);
    else
        f = pref.*(T1_ + T2_ + T3_ + T4_);
        g = f.'*Pe;
        G.yy(sel) = g(col);
    end
end
end

%% ------------------------------------------------------------------------
function R = stack_reflection(kzh, nh, kzs, ns, ds, is_te)
% Reflection coefficient of the stack ns (nearest layer first, last one a
% half-space) seen from the host layer, at the host/stack interface.
% s polarization (TM cylinders): r = (kz1 - kz2)/(kz1 + kz2)
% p polarization (TE cylinders, for H_y): r = (e2 kz1 - e1 kz2)/(e2 kz1 + e1 kz2)
nn = [nh, ns];
kk = [kzh, kzs];
R  = zeros(size(kzh));                       % nothing below the half-space
for j = numel(ns):-1:1                       % interface between layers j and j+1 of nn
    if is_te
        a = nn(j+1)^2*kk(:,j);  b = nn(j)^2*kk(:,j+1);
    else
        a = kk(:,j);            b = kk(:,j+1);
    end
    r = (a - b)./(a + b);
    if j == numel(ns)
        R = r;
    else
        ph = exp(2i*kk(:,j+1)*ds(j));        % two-way phase of the finite layer below this interface
        R  = (r + R.*ph)./(1 + r.*R.*ph);
    end
end
end

function s = follow_sqrt(arg, s_start)
% sqrt(arg) along a path, continuous from the value s_start before the first point
s = sqrt(arg);
flip = ones(size(s));
if abs(s(1) - s_start) > abs(s(1) + s_start), flip(1) = -1; end
j = abs(s(2:end) - s(1:end-1)) > abs(s(2:end) + s(1:end-1));
flip(find(j) + 1) = -1;
s = s.*cumprod(flip);
end

function [x, w] = panels(edges, xg, wg)
% Composite Gauss-Legendre rule on the panels between the edges
a = edges(1:end-1);  b = edges(2:end);
x = (a + b)/2 + (b - a)/2.*xg;               % nodes x panels
w = (b - a)/2.*wg;
x = x(:);  w = w(:);
end

function [x, w] = gauss_legendre(n)
% Nodes and weights on [-1, 1] (Golub-Welsch), as columns
b = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
w = 2*V(1, i).'.^2;
end
