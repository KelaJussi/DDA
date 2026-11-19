function ok = dda_layer_selftest()
% DDA_LAYER_SELFTEST  Checks of the layered background (dda_layer_green, dda_layer_field).
%
%   dda_layer_selftest      prints each check; returns true if all pass
%
% 1. A fictitious interface inside the host layer changes nothing.
% 2. Reciprocity: swapping source and observation point.
% 3. Layers with the host index give no reflected field.
% 4. A near-perfect mirror: the tangential field vanishes at the interface.
% 5. The field of a guided mode: far along a film the Green's function
%    advances with the effective index of the film's TE0 / TM0 mode.
% 6. The plane-wave field reproduces the Fresnel transmission coefficients.
% Checks 1-3 are also run at a complex frequency.
ok  = true;
lam = 0.6308e-6;  nh = 1.4665;
dx  = [0 0.2 1 5 20 80 -80 -5]*1e-6;
zi  = 0.10e-6*ones(size(dx));  zj = 0.15e-6*ones(size(dx));
rel = @(a, b) max(abs(a(:) - b(:)))/max(abs(b(:)));
pols = {'TM', 'TE'};

for k0 = 2*pi/lam*[1, 1 - 0.007i]
    tag = sprintf('(Im k0/Re k0 = %g)', imag(k0)/real(k0));
    La = struct('enable', true, 'd', 1.0e-6, 'nb', [1.33 1.52], 'db', 2e-6, 'na', [1.6 1], 'da', 0.3e-6);
    Lb = struct('enable', true, 'd', 0.6e-6, 'nb', [nh 1.33 1.52], 'db', [0.25e-6 2e-6], ...
                'na', [nh 1.6 1], 'da', [0.15e-6 0.3e-6]);
    L0 = struct('enable', true, 'd', 1.0e-6, 'nb', nh, 'db', [], 'na', nh, 'da', []);
    for p = 1:2
        Ga = dda_layer_green(k0, nh, La,  dx, zi + 0.25e-6, zj + 0.25e-6, pols{p});
        Gb = dda_layer_green(k0, nh, Lb,  dx, zi, zj, pols{p});
        Gc = dda_layer_green(k0, nh, La, -dx, zj + 0.25e-6, zi + 0.25e-6, pols{p});
        G0 = dda_layer_green(k0, nh, L0,  dx, zi, zj, pols{p});
        f = fieldnames(Ga);  e1 = 0;  e3 = 0;
        for n = 1:numel(f)
            e1 = max(e1, rel(Ga.(f{n}), Gb.(f{n})));
            e3 = max(e3, max(abs(G0.(f{n})(:))));
        end
        if p == 1, e2 = rel(Ga.yy, Gc.yy);
        else,      e2 = max([rel(Ga.xx, Gc.xx), rel(Ga.zz, Gc.zz), rel(Ga.xz, Gc.zx)]);
        end
        ok = report(ok, e1 < 1e-9, sprintf('%s fictitious interface %s', pols{p}, tag), e1);
        ok = report(ok, e2 < 1e-9, sprintf('%s reciprocity %s', pols{p}, tag), e2);
        ok = report(ok, e3 < 1e-12, sprintf('%s no interfaces %s', pols{p}, tag), e3);
    end
end

% 4. near-perfect mirror (index 60): direct + reflected tangential field at the interface
k0 = 2*pi/lam;  kh = nh*k0;
L  = struct('enable', true, 'd', 2e-6, 'nb', 60, 'db', [], 'na', nh, 'da', []);
x4 = [0.3 1]*1e-6;  z4 = 1e-10*[1 1];  s4 = 150e-9*[1 1];
R = sqrt(x4.^2 + (z4 - s4).^2);  u = kh*R;  cx = x4./R;
H0 = besselh(0,1,u);  H1 = besselh(1,1,u);  H2 = besselh(2,1,u);
G  = dda_layer_green(k0, nh, L, x4, z4, s4, 'TM');
e  = max(abs(1i/4*H0 + G.yy)./abs(1i/4*H0));
ok = report(ok, e < 0.05, 'TM mirror: E_y at the interface / direct field', e);
G  = dda_layer_green(k0, nh, L, x4, z4, s4, 'TE');
g0 = 1i/4*(H0 - H1./u + H2.*cx.*cx);
e  = max(abs(g0 + G.xx)./abs(g0));
ok = report(ok, e < 0.15, 'TE mirror: E_x at the interface / direct field', e);

% 5. guided modes of a film: air / 0.49 um film / 1.33, single mode per polarization
lam5 = 0.61309e-6;  k0 = 2*pi/lam5;  P = 0.435e-6;
L  = struct('enable', true, 'd', 0.49e-6, 'nb', 1.33, 'db', [], 'na', 1, 'da', []);
x5 = (150:183)*P;  z5 = 0*x5 + 1e-7;
G  = dda_layer_green(k0, nh, L, x5, z5, z5, 'TM');
ne = mean(angle(G.yy(2:end)./G.yy(1:end-1)) + 2*pi)/(k0*P);
n_te0 = slab_index(0.49e-6, lam5, nh, 1.33, 1, false);
ok = report(ok, abs(ne - n_te0) < 3e-3, sprintf('TE0 index of the film %.4f, from the Green''s function %.4f', n_te0, ne), abs(ne - n_te0));
G  = dda_layer_green(k0, nh, L, x5, z5, z5, 'TE');
ne = mean(angle(G.zz(2:end)./G.zz(1:end-1)) + 2*pi)/(k0*P);
n_tm0 = slab_index(0.49e-6, lam5, nh, 1.33, 1, true);
ok = report(ok, abs(ne - n_tm0) < 3e-3, sprintf('TM0 index of the film %.4f, from the Green''s function %.4f', n_tm0, ne), abs(ne - n_tm0));

% 6. Fresnel transmission, air -> host
k0 = 2*pi/lam;  kx = k0*sin(20*pi/180);
L  = struct('enable', true, 'd', 1e-6, 'nb', nh, 'db', [], 'na', 1, 'da', []);
k1 = sqrt(k0^2 - kx^2);  k2 = sqrt((nh*k0)^2 - kx^2);
F  = dda_layer_field(k0, nh, L, kx, 1e-6 - 1e-15, 'TM');
ok = report(ok, abs(abs(F.y) - 2*k1/(k1 + k2)) < 1e-6, 'plane wave, s polarization: Fresnel t', abs(abs(F.y) - 2*k1/(k1 + k2)));
F  = dda_layer_field(k0, nh, L, kx, 1e-6 - 1e-15, 'TE');
e  = abs(sqrt(abs(F.x)^2 + abs(F.z)^2) - 2*nh*k1/(nh^2*k1 + k2));
ok = report(ok, e < 1e-6, 'plane wave, p polarization: Fresnel t', e);

if ok, fprintf('All layer checks passed.\n'); else, fprintf('SOME LAYER CHECKS FAILED.\n'); end
end

function ok = report(ok, pass, name, value)
if pass, s = 'ok  '; else, s = 'FAIL'; ok = false; end
fprintf('%s  %-72s %.1e\n', s, name, value);
end

function n = slab_index(d, lam, nc, ns, na, tm)
% effective index of the fundamental mode of a slab (bisection on the dispersion relation)
k = 2*pi/lam;
    function v = f(n)
        kc = k*sqrt(nc^2 - n^2);  gs = k*sqrt(n^2 - ns^2);  ga = k*sqrt(n^2 - na^2);
        if tm, v = kc*d - atan(nc^2/ns^2*gs/kc) - atan(nc^2/na^2*ga/kc);
        else,  v = kc*d - atan(gs/kc) - atan(ga/kc);
        end
    end
a = max(ns, na) + 1e-9;  b = nc - 1e-9;
for it = 1:60
    n = (a + b)/2;
    if f(a)*f(n) <= 0, b = n; else, a = n; end
end
end
