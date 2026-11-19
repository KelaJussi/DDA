function ext = dda_extinction(P, Einc, Input, Array, subset)
% DDA_EXTINCTION  Extinction cross-section per unit area from dipole moments.
%
%   ext = dda_extinction(P, Einc, Input, Array)          whole array
%   ext = dda_extinction(P, Einc, Input, Array, subset)  only these particles
%
%   P, Einc  Ndof x Nk, components interleaved per particle (as in A)
%   subset   particle indices or logical mask (default: all)
%
%   ext(k) = 4*pi*k * sum Im(p . conj(E)) / (area * I0)
%     area = number of particles in the subset * unit-cell area
%     I0   = mean over the subset of |E|^2 summed over components
%   For SI polarizabilities (dda_alpha_is_SI: Mie and Lorentzian particles) the
%   factor 1/(4*pi*eps0*n^2) is applied, as in Calculate_3D_DDA_generic. Using I0 per particle makes 1D, 2D and 3D agree
%   (the old 3D code divided by 3 for the same reason; its 2D case did not).
%
% One formula for Calculate_DDA_eigenmodes, Calculate_3D_DDA_generic results
% (re-evaluated when saved) and the aperture plots.

N     = size(Array.r_i,1);
ncomp = size(P,1)/N;
Nk    = size(P,2);
if nargin < 5 || isempty(subset), subset = 1:N; end
if islogical(subset), subset = find(subset); end

Pp = reshape(double(P),    ncomp, N, Nk);
Ep = reshape(double(Einc), ncomp, N, Nk);
Pp = Pp(:,subset,:);
Ep = Ep(:,subset,:);

if any(Input.dimensionality == ["cyl_TM","cyl_TE"])
    cell_area = Input.p(1);                       % per unit length of the cylinders
else
    cell_area = Input.p(1)*Input.p(2);
end
area = numel(subset)*cell_area;
I0   = squeeze(mean(sum(abs(Ep).^2,1),2)).';    % 1 x Nk
work = squeeze(sum(imag(Pp.*conj(Ep)),[1 2])).'; % 1 x Nk

k   = double(Array.k(:).');
ext = 4*pi*k.*work./(area*I0);

if dda_alpha_is_SI(Input, Array)
    eps0 = 8.8541878188e-12;
    n    = double(Array.n(:).');
    ext  = ext/4/pi/eps0./n.^2;
end
end
