function [A]=Make_1D_int_matrix_A_sub(Input,Array,kk)
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

gain = Input.gain;
exp_term=exp(1i*k.*R)./(R);   
%exp_term=exp(1i*k.*R)./(4*pi*R);   
R_sqrd_term=(1i*k.*R-1)./(R.^2);
%gain_term = exp(gain.*k.*R);

r_i_per_R(:,:,1)=(r_i(:,1)-r_i(:,1).')./R;
r_i_per_R(:,:,2)=(r_i(:,2)-r_i(:,2).')./R;
r_i_per_R(:,:,3)=(r_i(:,3)-r_i(:,3).')./R;
%% Calculate off-diagonal blocks for matrix Aij for single wavelength


switch Input.dimensionality
    case '1Dx'  
        ind_ii=1;
        ind_jj=1;
        
        I3=eye(3);
    case '1Dy'
        ind_ii=2;
        ind_jj=2;
        
        I3=eye(3);
    case '1Dz'
        ind_ii=3;
        ind_jj=3;
         I3=eye(3);
end

% Form the A and Aself matrices:
% A=zeros(size(R_fun(Array.r_i),1),size(R_fun(Array.r_i),1));
eps0 =  8.8541878188e-12; %F/m
if 0%(ind_ii==3 && ind_jj==3)
    %% A only function of R etc.
    A = exp(1i*k.*R_fun(Array.r_i))./(4*pi*R_fun(Array.r_i)).*(k^2.*(-1)+...
        (1i*k.*R_fun(Array.r_i)-1)./(R_fun(Array.r_i).^2).*(-1))+...
        diag(repmat(squeeze(1./Array.alpha(1,3,3,kk)),Input.N,1));       % A_ii , i.e. polarizabilities
    Aself=  diag(4*pi*eps0./Array.alpha(:,ind_ii,ind_jj,kk));
else
    %% R_hats are needed in calculation of A
    for ii=ind_ii     % loop over the right 3 x 3 indices of the dyadic
        for jj=ind_jj
            % A = exp(1i*k.*R_fun(Array.r_i))./(4*pi*R_fun(Array.r_i)).* ...
            %     (k^2.*(R_jk_hat_fun(Array.r_i,ii).*R_jk_hat_fun(Array.r_i,jj)-I3(ii,jj)) ...
            %     +(1i*k.*R_fun(Array.r_i)-1)./(R_fun(Array.r_i).^2).*...
            %     (3*R_jk_hat_fun(Array.r_i,ii).*R_jk_hat_fun(Array.r_i,jj))-I3(ii,jj))+...
            %     diag(repmat(squeeze(1./Array.alpha(1,ind_ii,ind_jj,kk)),Input.N,1));
             % AB = exp(1i*k.*R_fun(Array.r_i))./(R_fun(Array.r_i)).* ...
             %    (k^2.*(R_jk_hat_fun(Array.r_i,ii).*R_jk_hat_fun(Array.r_i,jj)-I3(ii,jj)) ...
             %    +(1i*k.*R_fun(Array.r_i)-1)./(R_fun(Array.r_i).^2).*...
             %    (3*R_jk_hat_fun(Array.r_i,ii).*R_jk_hat_fun(Array.r_i,jj))-I3(ii,jj));%+...

             A = exp_term.*(k^2.*(r_i_per_R(:,:,ii).*r_i_per_R(:,:,jj)-I3(ii,jj)) ...
                  +R_sqrd_term.*(3.*r_i_per_R(:,:,ii).*r_i_per_R(:,:,jj)-I3(ii,jj))); 
            %ABC = A-AB;
             % Aself=  diag(4*pi*eps0./Array.alpha(:,ind_ii,ind_jj,kk));
              Aself=  diag(1./double(Array.alpha(:,ind_ii,ind_jj,kk)));
        end
    end
end
nnzNaN=find(isnan(A));
n=1;
%A=A+Aself;
if dda_alpha_is_SI(Input,Array)   % SI polarizabilities: Gaussian-form terms / (4*pi*eps0*n^2)
        A = A/4/pi./eps0./Array.n(kk).^2;
end
A((nnzNaN))=1e-10;
%A((nnzNaN))=0;
%A((nnzNaN))=0;
A(isinf(A))=1e-10;  % 1/1e-30
%A(isinf(A))=0;  % 1/1e-30
if Input.particle.r == 0
    Aself = zeros(size(Aself));
elseif Input.particle.r == 1
    Aself = 1e29*diag(ones(size(Aself,1),1));
end
A=A+Aself;


end