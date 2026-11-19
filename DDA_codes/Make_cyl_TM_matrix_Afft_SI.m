function [Afft,siz]=Make_cyl_TM_matrix_Afft_SI(kk,Array, r_i, Input)
% Form A = N x N block interaction matrix, where i=j terms are inverses of
% single particle polarizabilities alpha

%r_i should now be r_i_conv marking all the possible distances between
%particles in evenly periodic array.
    
kh = Array.k(kk);
%% 
%r_i = Array.r_i;
%R = sqrt((r_i(:,1)-r_i(:,1).').^2 + (r_i(:,2)-r_i(:,2).').^2 + (r_i(:,3)-r_i(:,3).').^2); %With cylinders, y shouldn't have particles (only 1) so only 2 possible dimensions ,
R=sqrt(r_i(:,1).^2 + r_i(:,2).^2 + r_i(:,3).^2);
% although only one is used in practice with cylinders. Still 3d is used
% for consistency

R(R<1e-14)=0;

%% Calculate off-diagonal blocks for matrix Aij for single wavelength


% Form the A and Aself matrices:
% A=zeros(size(R_fun(Array.r_i),1),size(R_fun(Array.r_i),1));
eps0 =  8.8541878188e-12; %F/m

%% Cylinder TE scalar 
N = Input.N;
A = zeros(N,N,'like',Array.r_i);

%R = R_fun(Array.r_i);   % NxN distances
% 2D Green function
G = 1i/4 * besselh(0,1,kh*R);

% remove self-interaction for now
G(1:N+1:end) = 0;

% Fill interaction matrix
A = -kh.^2 * G ./ (eps0.*Array.n(kk).^2);
nnzNaN=find(isnan(A));

%A=A+Aself;

A((nnzNaN))=1e-10;

A(isinf(A))=1e-10;  % 1/1e-30

A_BT=reshape(A,[Input.Nmax(1)*2-1 Input.Nmax(2)*2-1  Input.Nmax(3)*2-1]);
          

%% Check the sizes for making the correct ffts
if Input.Nmax(3)> 1
    siz.A = size(A_BT);
else
    siz.A = [size(A_BT) 1];
end
siz.X = [Input.Nmax(1) Input.Nmax(2) Input.Nmax(3)];
siz.total =  2.^nextpow2(siz.A + siz.X -1);
% siz.total =  siz.A + siz.X -1;

%% Fourier and inverse fourier transforms to perform the convolution:
Afft=cell(1,1);
Afft=(fftn(A_BT,siz.total));

% Add self term (polarizability)
%for i = 1:N
%    A(i,i) = 1 / Array.alpha(i,kk);   % scalar alpha now
%end

% Extract unique distances (exclude zeros)
%{
Rvec = R(:);
Gvec = G(:);

mask = Rvec > 1e-12;   % remove self terms
Rvec = Rvec(mask);
Gvec = Gvec(mask);

% Sort by distance
[Rvec_sorted, idx] = sort(Rvec);
Gvec_sorted = Gvec(idx);

figure;
subplot(2,1,1)
plot(Rvec_sorted, real(Gvec_sorted), 'b', 'LineWidth', 1.5); hold on
plot(Rvec_sorted, imag(Gvec_sorted), 'r', 'LineWidth', 1.5);
xlabel('Distance R')
ylabel('G(R)')
legend('Re(G)', 'Im(G)')
title('2D Green function G = i/4 H_0^{(1)}(k_h R)')
grid on

subplot(2,1,2)
plot(Rvec_sorted, abs(Gvec_sorted), 'k', 'LineWidth', 1.5);
xlabel('Distance R')
ylabel('|G(R)|')
title('|G(R)|')
grid on
%}

% A = -k.^2./(eps0.*Array.n(kk).^2).*...
%     exp(1i*k*R)./(4*pi*R) .* ((1 + 1i./(k*R) - 1./(k*R).^2)*I3 ... 
 %    +(-1 - 3i./(k*R) + 3./(k*R).^2) .* RhatRhat);

  %Aself=  diag(1./double(Array.alpha(:,ind_ii,ind_jj,kk)));





end