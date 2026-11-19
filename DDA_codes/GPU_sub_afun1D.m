function [y]=GPU_sub_afun1D(x,Afft,siz,alpha,Input)
%% Performs fast matrix-vector multiplication related to iteratively solving
% source polarization in DDA:  
%             Ax, where A = interaction matrix 
%                       x = trial of source polarization
%
% Script is based on the fact that A is Block Toeplitz (BT) for rectangular
% dipole arrays, and hence the multiplication can be turned into a
% convolution. The convolution is calculated using the regular fft-based
% speed-up schemes. More info can be found from e.g. 
%
% Draine papers, Flatau Review on DDA 
%
% In addition, only small part of the A has to computed and
% convolved as there's a lot of redundancy in the full A matrix for the
% system. 
%
% correct size for the interaction matrix is 2N-1, if the traditional A is NxN
% But we only need to calculate a part of A, where the dipole has N-1 neighbours in 
% all directions. 
% e.g. N=10, and NxN array traditionally requires constructing A of size N^2xN^2 
%
% Things to make into the codes, 
% 1) check that indeed we're working with square/rectangular arrays
% Input.arrayN  or Input.N
%
% How the actual algorithm works:
% 1) construct the right A_block for convolution
% 2) check sizes and make ffts
% 3) make ifft and crop the right part of y
%
% 1) For NxN array, we need to calculate A_block = B, which corresponds to 
%    an interaction matrix A of a dipole surrounded by N-1 neighbouring dipoles
%    in all directions.
%
% x=gpuArray(single(x));
if nnz(strfind(Input.inversion_tool,'gpu'))>0
    x=gpuArray(single(x));
end

switch Input.dimensionality
    case '1Dx'  
        Afft = Afft{1,1};                        
    case '1Dy'
        Afft = Afft{2,2};    
    case '1Dz'
        Afft = Afft{3,3};  
    otherwise
        ind_ii=1:3;      ind_jj=1:3; 
end

%% 2.  structure x:
% Zero padding and reshaping of x in to X 
Ntot=Input.Nmax(1)*Input.Nmax(2)*Input.Nmax(3); 
X_BT=single(reshape(x(1:Ntot),Input.Nmax(1),Input.Nmax(2),Input.Nmax(3)));

%X_BT=single(reshape(x(1:Ntot),Input.Nmax(2),Input.Nmax(1),Input.Nmax(3)));


 
%% 3. Fourier and inverse fourier transforms to perform the convolution:
% OBS! fftn might be able to do an even higher dimensionality FFT faster than with now is being for looped! 10.12.2015 mjh
y=zeros(size(x),'like',alpha);
%y_full=cell(length(ind_jj),1);

Xfft=fftn(X_BT,siz.total);  % (jj,x)

%% 4. Calculate matrix Aij for single wavelength                  
%y_full=ifftn(Afft{ii,jj}.*fftn(X_BT{ii},siz.total));          % actual convolution part
y_full=ifftn(Afft.*Xfft);          % actual convolution part
y=reshape(y_full(ceil((siz.A(1)-1)/2)+(1:siz.X(1)),ceil((siz.A(2)-1)/2)+(1:siz.X(2)),:),[],1);
 %y=reshape(y_full(ceil((siz.A(2)-1)/2)+(1:siz.X(2)),ceil((siz.A(1)-1)/2)+(1:siz.X(1)),:),[],1);
 

%% 5. Form y = alpha_i^-1 * P_i   +   A_ij * P_j , which is the full trial
y = (1./alpha).*x + y;

end
