function [A]=Make_int_matrix_A_SI(Input,Array,kk)
% Form A = N x N block interaction matrix, where i=j terms are inverses of
% single particle polarizabilities alpha
% Made by Claude, Original from JTK in folder as Make_int_matrix_A_SI - JTK Original
%
% Optional boundary settings:
%   Input.mirror_A          [s1 s2 s3 s4] image weights for mirrors at
%                           x = 0, x = max(x), y = 0, y = max(y)  (0 = off)
%   Input.periodic_boundary [px py] logical, periodic boundary in x and/or y.
%                           Rectangular cell, no Bloch phase (k_parallel = 0).
%   Input.periodic_cell     [Lx Ly] unit cell size in metres (the GUI field is in nm).
%                           Optional: if missing, it is inferred from the
%                           particle positions as extent + nearest-neighbour
%                           spacing (assumes a uniformly spaced grid). It must
%                           be larger than the array's extent in that direction.
%   Input.periodic_Ncells   [Nx Ny] number of image cells summed on EACH side
%                           of the central cell (default [10 10]): 0 = no
%                           images, 1 = the cells at -L and +L (3 cells in all),
%                           10 = 21 cells.
%   A periodic direction cannot be combined with a mirror in the same direction.


k=Array.k(kk);
%% JTK help funcitons
r_i = double(Array.r_i);   % positions are stored single: the matrix must be built in double
R = sqrt((r_i(:,1)-r_i(:,1).').^2 + (r_i(:,2)-r_i(:,2).').^2 + (r_i(:,3)-r_i(:,3).').^2);
R(R<1e-14)=0;
eps0 =  8.8541878188e-12; %F/m

%% Boundary options and consistency checks
mirrorA = zeros(1,4);
if isfield(Input,'mirror_A')
    mirrorA(1:numel(Input.mirror_A)) = Input.mirror_A;
end

pbc = [false false];
if isfield(Input,'periodic_boundary')
    if numel(Input.periodic_boundary) ~= 2
        error('Make_int_matrix_A_SI:periodicSize', ...
            'Input.periodic_boundary must have two elements [x, y].');
    end
    pbc = logical(Input.periodic_boundary(:).');
end

if pbc(1) && any(mirrorA(1:2))
    error('Make_int_matrix_A_SI:periodicMirrorX', ...
        ['Periodic boundary in x (Input.periodic_boundary(1)) cannot be used together ' ...
         'with an x mirror (Input.mirror_A(1:2) = [%g %g]).'], mirrorA(1), mirrorA(2));
end
if pbc(2) && any(mirrorA(3:4))
    error('Make_int_matrix_A_SI:periodicMirrorY', ...
        ['Periodic boundary in y (Input.periodic_boundary(2)) cannot be used together ' ...
         'with a y mirror (Input.mirror_A(3:4) = [%g %g]).'], mirrorA(3), mirrorA(4));
end

%% Periodic lattice translations
% Tall includes the central cell (0,0); Tper excludes it.
if any(pbc)
    Lcell = get_periodic_cell(Input, r_i, pbc);
    Ncells = [10 10];
    if isfield(Input,'periodic_Ncells')
        Ncells = double(Input.periodic_Ncells(:).') .* [1 1]; % scalar -> [N N]
    end
    nx = 0; ny = 0;
    if pbc(1), nx = -Ncells(1):Ncells(1); end
    if pbc(2), ny = -Ncells(2):Ncells(2); end
    [NX,NY] = ndgrid(nx,ny);
    Tall = [NX(:)*Lcell(1), NY(:)*Lcell(2), zeros(numel(NX),1)];
    Tper = Tall(NX(:)~=0 | NY(:)~=0, :);
    overlap_tol = 1e-3*min(Input.p(1:2));
else
    Tall = [0 0 0];
end

%% Calculate off-diagonal blocks for matrix Aij for single wavelength
Rhat(:,:,1) = ((r_i(:,1)-r_i(:,1).')./R);
Rhat(:,:,2) = ((r_i(:,2)-r_i(:,2).')./R);
 Rhat(:,:,3) = ((r_i(:,3)-r_i(:,3).')./R);
switch Input.dimensionality
    case '1Dx'  
        ind = 1;
        ind_ii=1;
        ind_jj=1;
        A = zeros(size(r_i,1),size(r_i,1),'like',r_i); % y=zeros(size(x),'like',alpha);
 
    case '1Dy'
        ind = 1;
        ind_ii=2;
        ind_jj=2; 
        A = zeros(size(r_i,1),size(r_i,1),'like',r_i); % y=zeros(size(x),'like',alpha);
 
    case '1Dz'
        ind = 1;
        ind_ii=3;
        ind_jj=3;
          A = zeros(size(r_i,1),size(r_i,1),'like',r_i); % y=zeros(size(x),'like',alpha);
 
    case '2D'
        ind = 2;
        ind_ii = 1:ind;
        ind_jj = 1:ind;
        A = zeros(ind*size(r_i,1),ind*size(r_i,1),'like',r_i); % y=zeros(size(x),'like',alpha);
    case '3D'
        ind = 3;
        ind_ii = 1:ind;
        ind_jj = 1:ind;
        A = zeros(ind*size(r_i,1),ind*size(r_i,1),'like',r_i); % y=zeros(size(x),'like',alpha);
        
        
end
Aself = A;
I3=eye(3);
off = ind_ii(1) - 1; % block index offset: 0 for 2D/3D, needed for '1Dy'/'1Dz'

for ii = ind_ii
    for jj = ind_jj
         A (ii-off:ind:end,jj-off:ind:end)= -k.^2./(eps0.*Array.n(kk).^2).*...
             exp(1i*k*R)./(4*pi*R) .* ((1 + 1i./(k*R) - 1./(k*R).^2)*I3(ii,jj) ... 
             +(-1 - 3i./(k*R) + 3./(k*R).^2) .* Rhat(:,:,ii) .* Rhat(:,:,jj));

    end
end

if isfield(Array, 'alpha') || isfield(Array, 'alpha_types')
    % Self term = inverse of each particle's polarizability tensor. This must be
    % a matrix inverse: an element-wise 1./alpha turns the near-zero
    % off-diagonal entries (1e-30 placeholders in Mie_polarizability_DDA) into
    % couplings ~1e40 that swamp the physics in 2D and 3D. (Same as MJH's
    % Make_3D_int_matrix_A_sub, which uses alpha\eye(3).)
    % Works with full and compact storage (dda_alpha): with compact storage each
    % particle type is inverted once.
    if isfield(Array, 'alpha_types')
        types = Array.alpha_type;  M = size(Array.alpha_types, 1);
        Ainv = cell(1, M);
        for m = 1:M
            blk = double(squeeze(Array.alpha_types(m,ind_ii,ind_jj,kk)));
            if ind == 1, Ainv{m} = 1./blk; else, Ainv{m} = inv(blk); end
        end
    end
    for ii = 1:Input.N
        rows = (ind*ii-ind+1):(ind*ii);
        if isfield(Array, 'alpha_types')
            Aself(rows,rows) = Ainv{types(ii)};
        else
            blk = double(squeeze(Array.alpha(ii,ind_ii,ind_jj,kk)));
            if ind == 1, Aself(rows,rows) = 1./blk; else, Aself(rows,rows) = inv(blk); end
        end
    end
    Aself(isinf(Aself))=1e-10;
else
    Aself = 0;
end

nnzNaN=find(isnan(A));

A((nnzNaN))=1e-10;

A(isinf(A))=1e-10;  % 1/1e-30

nullresponse = -1i*k.^3./(6*pi*eps0*Array.n(kk).^2);
Aself(1,1);
A=A+Aself;

%% Periodic images of the array (central cell is already in A)
if any(pbc)
    A = A + green_lattice_sum(k, Array.n(kk), r_i, r_i, Tper, ind, ind_ii, ind_jj, overlap_tol);
end

%% Mirror images
% If the other direction is periodic, each mirror image is also repeated
% along that direction (Tall). Without periodicity Tall = [0 0 0].
% The image of a dipole p is s*M*p, with M the reflection (the component
% normal to the mirror changes sign) and s = mirror_A:
%   2D/3D  s = +1: field mirror-symmetric (tangential components even, normal
%                  odd; a magnetic-wall / PMC plane)
%          s = -1: field mirror-antisymmetric (tangential odd, normal even; a
%                  perfectly conducting / PEC plane)
%   1Dx/1Dy/1Dz (one component): the scalar meaning p(-x) = s*p(x), as before.
% The image coupling is therefore s * G(r_i - M r_j) * M.
%
% The planes go through the outermost particles (x = 0, x = max(x), ...), so
% the particles ON a plane are their own images. Their images are not added:
% that counted the coupling to them twice (three times for a single column
% between two mirrors), which made A non-passive: negative extinction
% (transmission > 1) and quasi-normal modes with Im(omega) > 0.
% The components that the symmetry makes odd on the plane (s*M = -1) must be
% zero there; they are decoupled from the other particles (p = alpha*Einc,
% zero for an incident field with the same symmetry).
if ind > 1
    Mrefl = {[-1 1 1], [1 -1 1]};         % reflection of the dipole components (x, y mirror)
else
    Mrefl = {[1 1 1], [1 1 1]};           % one component: scalar sign only
end
planes  = [0, max(r_i(:,1)), 0, max(r_i(:,2))];   % [x=0, x=max(x), y=0, y=max(y)]
on_tol  = 1e-3*min(Input.p(1:2));
decouple = false(size(A,1), 1);
for m = find(mirrorA)
    d = 1 + (m > 2);                      % 1: x mirror, 2: y mirror
    M = Mrefl{d};
    r_img = r_i;
    r_img(:,d) = 2*planes(m) - r_i(:,d);
    Gm = green_lattice_sum(k, Array.n(kk), r_i, r_img, Tall, ind, ind_ii, ind_jj, [], M);
    on = abs(r_i(:,d) - planes(m)) < on_tol;
    if any(on)
        dofs = particle_dofs(find(on), ind);          %#ok<FNDSB>  ind x n_on
        Gm(:, dofs(:)) = 0;                           % no second copy of on-plane particles
        odd = mirrorA(m)*M(ind_ii) < 0;               % components that must vanish on the plane
        decouple(reshape(dofs(odd,:), [], 1)) = true;
    end
    A = A + mirrorA(m)*Gm;
end
if any(decouple)
    dA = diag(A);
    A(decouple, :) = 0;
    A(:, decouple) = 0;
    A(sub2ind(size(A), find(decouple), find(decouple))) = dA(decouple);
end
%% ELA but actually they are big particles (perfectly polarizable but with varying real parts)
%{

if Input.particle.r == 0
    %Aself = zeros(size(Aself));
    nullresponse = -1i*k.^3./(6*pi*eps0*Array.n(kk).^2);
     Aself = nullresponse.*diag(ones(size(Aself,1),1));
elseif Input.particle.r == 1e-9
    nullresponse = -1i*k.^3./(6*pi*eps0*Array.n(kk).^2);
    Aself = (1e29+nullresponse)*diag(ones(size(Aself,1),1));
elseif Input.particle.r == 2e-9
    nullresponse = -1i*k.^3./(6*pi*eps0*Array.n(kk).^2);
    Aself = (1e30+nullresponse)*diag(ones(size(Aself,1),1));
elseif Input.particle.r == 3e-9
    nullresponse = -1i*k.^3./(6*pi*eps0*Array.n(kk).^2);
    Aself = (5e29+nullresponse)*diag(ones(size(Aself,1),1));
elseif Input.particle.r == 4e-9
    nullresponse = -1i*k.^3./(6*pi*eps0*Array.n(kk).^2);
    Aself = (8e29+nullresponse)*diag(ones(size(Aself,1),1));
elseif Input.particle.r == 5e-9
    nullresponse = -1i*k.^3./(6*pi*eps0*Array.n(kk).^2);
    Aself = (7e29+nullresponse)*diag(ones(size(Aself,1),1));
elseif Input.particle.r == 6e-9
    nullresponse = -1i*k.^3./(6*pi*eps0*Array.n(kk).^2);
    Aself = (3e29+nullresponse)*diag(ones(size(Aself,1),1));
elseif Input.particle.r == 7e-9
    nullresponse = -1i*k.^3./(6*pi*eps0*Array.n(kk).^2);
    Aself = (2e29+nullresponse)*diag(ones(size(Aself,1),1));
end
%}





end

%% ------------------------------------------------------------------------
function dofs = particle_dofs(idx, ind)
% Matrix rows of particles idx: ind x numel(idx), components interleaved
dofs = (idx(:).' - 1)*ind + (1:ind).';
end

%% ------------------------------------------------------------------------
function G = green_lattice_sum(k, n_med, r_obs, r_src, T, ind, ind_ii, ind_jj, overlap_tol, colsign)
% Sum of the same dyadic Green's function blocks as in the main function,
% G(i,j) = sum_t G0(r_obs(i,:) - r_src(j,:) - T(t,:)) * diag(colsign).
% colsign (default [1 1 1]) multiplies the source dipole components: the
% reflection M for mirror images.
% NaN/Inf entries (coinciding points, e.g. particles on a mirror plane) are
% set to 1e-10 as in the original code. If overlap_tol is non-empty, a
% coinciding point raises an error instead (used for periodic images).
eps0 = 8.8541878188e-12; %F/m
if nargin < 10 || isempty(colsign), colsign = [1 1 1]; end
I3   = eye(3);
off  = ind_ii(1) - 1;
G    = zeros(ind*size(r_obs,1), ind*size(r_src,1), 'like', r_obs);
pref = -k.^2./(eps0.*n_med.^2);

for t = 1:size(T,1)
    RX = r_obs(:,1) - (r_src(:,1).' + T(t,1));
    RY = r_obs(:,2) - (r_src(:,2).' + T(t,2));
    RZ = r_obs(:,3) - (r_src(:,3).' + T(t,3));
    Rt = sqrt(RX.^2 + RY.^2 + RZ.^2);

    if ~isempty(overlap_tol) && any(Rt(:) < overlap_tol)
        error('Make_int_matrix_A_SI:periodicOverlap', ...
            ['A particle coincides with a periodic image of a particle ' ...
             '(cell translation [%g %g]). Check Input.periodic_cell.'], T(t,1), T(t,2));
    end

    kR = k*Rt;
    g  = pref .* exp(1i*kR)./(4*pi*Rt);
    f1 = 1 + 1i./kR - 1./kR.^2;
    f2 = -1 - 3i./kR + 3./kR.^2;
    Rh = {RX./Rt, RY./Rt, RZ./Rt};

    for ii = ind_ii
        for jj = ind_jj
            blk = g .* (f1*I3(ii,jj) + f2 .* Rh{ii} .* Rh{jj}) * colsign(jj);
            blk(isnan(blk) | isinf(blk)) = 1e-10;
            G(ii-off:ind:end, jj-off:ind:end) = G(ii-off:ind:end, jj-off:ind:end) + blk;
        end
    end
end
end

%% ------------------------------------------------------------------------
function Lcell = get_periodic_cell(Input, r_i, pbc)
% Unit cell size [Lx Ly] (zero in non-periodic directions).
Lcell = [0 0];
if isfield(Input,'periodic_cell') && ~isempty(Input.periodic_cell)
    Lcell = double(Input.periodic_cell(:).') .* [1 1]; % scalar -> [L L]
end
Lcell(~pbc) = 0;
dirs = 'xy';
for d = find(pbc)
    if Lcell(d) <= 0
        % Infer: particle extent + nearest-neighbour spacing (uniform grid)
        x    = sort(r_i(:,d));
        span = x(end) - x(1);
        dx   = diff(x);
        dx   = dx(dx > 1e-6*span);
        if span == 0 || isempty(dx)
            error('Make_int_matrix_A_SI:periodicCell', ...
                ['Cannot infer the periodic cell size in %s from the particle ' ...
                 'positions. Set Input.periodic_cell = [Lx Ly].'], dirs(d));
        end
        Lcell(d) = span + min(dx);
    end
    % The images of the array must not overlap it: e.g. 700e-9 typed into the
    % GUI's nm field gives 7e-16 m, which put every image on top of a particle.
    x    = r_i(:,d);
    span = max(x) - min(x);
    tol  = 1e-3*min(Input.p(1:2));
    if Lcell(d) <= span + tol
        error('Make_int_matrix_A_SI:periodicCellSize', ...
            ['Periodic cell size in %s is %.3g m, but the array is %.3g m long in %s: the ' ...
             'periodic images would overlap the particles. Input.periodic_cell is in metres ' ...
             '(the GUI field is in nm), e.g. [700e-9 0] in code or 700 in the GUI.'], ...
            dirs(d), Lcell(d), span, dirs(d));
    end
end
end