function [A]=Make_3D_int_matrix_A_sub_JTK(Input,Array,kk)
% Form A = N x N block interaction matrix, where i=j terms are inverses of
% single particle polarizabilities alpha
    
%% For calculating R
    function R=R_fun(r_i)
        R=sqrt(Array.r_ij_fun(r_i,1).^2 + Array.r_ij_fun(r_i,2).^2);   % r_z == 0 identically in XY-plane
        
        % replace zero R terms with large distance so, that those won't affect results             
        R(1:size(R,1)+1:end) = inf;     
    end
%% For calculating R_hat
    function R_jk_hat=R_jk_hat_fun(r_i,ii)
        if ii<3
            R_jk_hat =  Array.r_ij_fun(r_i,ii)./R_fun(r_i);
        else
            R_jk_hat=zeros(size(r_i,1));      % r_z == 0 identically in XY-plane
        end
    end
k=Array.k(kk);
%% JTK help funcitons
r_i = Array.r_i;
R=sqrt((r_i(:,1)-r_i(:,1).').^2 + (r_i(:,2)-r_i(:,2).').^2 + (r_i(:,3)-r_i(:,3).').^2);
R(R<1e-14)=0;

exp_term=exp(1i*k.*R)./(R);   
%exp_term=exp(1i*k.*R)./(4*pi*R);   
R_sqrd_term=(1i*k.*R-1)./(R.^2);

r_i_per_R(:,:,1)=(r_i(:,1)-r_i(:,1).')./R;
r_i_per_R(:,:,2)=(r_i(:,2)-r_i(:,2).')./R;
r_i_per_R(:,:,3)=(r_i(:,3)-r_i(:,3).')./R;
%% Calculate off-diagonal blocks for matrix Aij for single wavelength


ind_ii=1:3;
ind_jj = 1:3;
I3=eye(3);

% Form the A and Aself matrices:
% A=zeros(size(R_fun(Array.r_i),1),size(R_fun(Array.r_i),1));
eps0 =  8.8541878188e-12; %F/m
A=zeros(3*size(R_fun(Array.r_i),1),3*size(R_fun(Array.r_i),1),'like',Array.r_i); % y=zeros(size(x),'like',alpha);
    
%% R_hats are needed in calculation of A
for ii=ind_ii     % loop over the right 3 x 3 indices of the dyadic
    for jj=ind_jj


         A(ii:3:end,jj:3:end) = exp_term.*(k^2.*(r_i_per_R(:,:,ii).*r_i_per_R(:,:,jj)-I3(ii,jj)) ...
              +R_sqrd_term.*(3.*r_i_per_R(:,:,ii).*r_i_per_R(:,:,jj)-I3(ii,jj))); 

    end
end
for ii=1:Input.N  
        Aself((3*ii-2):(3*ii),(3*ii-2):(3*ii))=squeeze(Array.alpha(ii,:,:,kk))\eye(3);   % matrix inverse, not 1./alpha
    %     Aself((3*ii-2):(3*ii),(3*ii-2):(3*ii))=tmp;
end

%A=A+Aself;
if dda_alpha_is_SI(Input,Array)   % SI polarizabilities: Gaussian-form terms / (4*pi*eps0*n^2)
        A = A/4/pi./eps0./Array.n(kk).^2;
end
[J]=find(Aself);    % find non-zero self interaction terms and input into the A(ij)
    A(J)=Aself(J);
nnzNaN=find(isnan(A));

A((nnzNaN))=1e-10;
%A((nnzNaN))=0;
A(isinf(A))=0;  %1e-10 ,1/1e-30
%A=A+Aself;


end