function [A,Aself]=Make_3D_int_matrix_A_sub(Input,Array,kk)
    % Form A = 3N x 3N block interaction matrix, where i=j terms are inverses of
    % single particle polarizabilities alpha
    
    % %% make r_ij and R functions for A calculations
    % Array.r_ij_fun = @(r_i,ii) repmat(r_i(:,ii),[1 size(r_i,1)])-...  % x-component
    %                            repmat(r_i(:,ii).',[size(r_i,1) 1]);   % y-component
    %        
    % Array.R_fun = @(r_i) sqrt(Array.r_ij_fun(r_i,1).^2 + Array.r_ij_fun(r_i,2).^2);   % r_z == 0 identically in XY-plane
    
    %% nested functions to calculate the R and hat_R
        function R=R_fun(r_i)
            R=sqrt(Array.r_ij_fun(r_i,1).^2 + Array.r_ij_fun(r_i,2).^2);   % r_z == 0 identically in XY-plane
            
            I=find(R==0);   % replace zero R terms with large distance so, that those won't affect results
            R(I)=1e6;        
        end
    
        function R_jk_hat=R_jk_hat_fun(r_i,ii)
            if ii<3
    %             R_jk_hat{ii} = Array.r_ij{ii}./R;
                R_jk_hat =  Array.r_ij_fun(r_i,ii)./R_fun(r_i);
            else
    %             R_jk_hat{3} = zeros(size(R_jk_hat{1}));  
                R_jk_hat=zeros(size(r_i,1));                                   % r_z == 0 identically in XY-plane
            end
        end
    
    k=Array.k(kk);
         
    %% Calculate off-diagonal blocks for matrix Aij for single wavelength
    
    % Form the A and Aself matrices:
    A=zeros(3*size(R_fun(Array.r_i),1),3*size(R_fun(Array.r_i),1),'like',Array.r_i); % y=zeros(size(x),'like',alpha);
    
    Aself=A;
    
    I3=eye(3);
    if strfind(Input.inversion_tool,'gpu')
        I3=gpuArray(single(I3));
    end
    k=Array.k(kk);
    
    % if kk==1
    %     save tmp_gpu2
    % end
    % tic
    R=sqrt(Array.r_ij_fun(Array.r_i,1).^2 + Array.r_ij_fun(Array.r_i,2).^2);
    R_jk_hat{1} =  Array.r_ij_fun(Array.r_i,1)./R;
    R_jk_hat{2} =  Array.r_ij_fun(Array.r_i,2)./R;
    R_jk_hat{3} =  zeros(size(Array.r_i,1));          
    % Array.r_ij_fun = @(r_i,ii) repmat(r_i(:,ii),[1 size(r_i,1)])-...  % x-component
    %                            repmat(r_i(:,ii).',[size(r_i,1) 1]);   % y-component
    %        
    % Array.R_fun = @(r_i) sqrt(Array.r_ij_fun(r_i,1).^2 + Array.r_ij_fun(r_i,2).^2);   % r_z == 0 identically in XY-plane
    
    for ii=1:3     % loop over the 3 x 3 indices of the dyadic 
        for jj=1:3        
            A(ii:3:end,jj:3:end) = exp(1i*k.*R)./(R).* ...
                                   (k^2.*(R_jk_hat{ii}.*R_jk_hat{jj}-I3(ii,jj)) ...
                                     +(1i*k.*R-1)./(R.^2).*...
                                      (3*R_jk_hat{ii}.*R_jk_hat{jj})-I3(ii,jj));
           
    %         A(ii:3:end,jj:3:end) = exp(1i*k.*R_fun(Array.r_i))./(R_fun(Array.r_i)).* ...
    %                                (k^2.*(R_jk_hat_fun(Array.r_i,ii).*R_jk_hat_fun(Array.r_i,jj)-I3(ii,jj)) ...
    %                                  +(1i*k.*R_fun(Array.r_i)-1)./(R_fun(Array.r_i).^2).*...
    %                                   (3*R_jk_hat_fun(Array.r_i,ii).*R_jk_hat_fun(Array.r_i,jj))-I3(ii,jj));
                                  
    %         A(ii:3:end,jj:3:end) = exp(1i*k.*R_fun(Array.r_i))./(4*pi*R_fun(Array.r_i)).* ...
    %                                (k^2.*(R_jk_hat_fun(Array.r_i,ii).*R_jk_hat_fun(Array.r_i,jj)-I3(ii,jj)) ...
    %                                  +(1i*k.*R_fun(Array.r_i)-1)./(R_fun(Array.r_i).^2).*...
    %                                   (3*R_jk_hat_fun(Array.r_i,ii).*R_jk_hat_fun(Array.r_i,jj))-I3(ii,jj));                      
        end
    end
    if dda_alpha_is_SI(Input,Array)   % SI polarizabilities (was: Au_Mie/glass_Mie only, and without n^2)
        eps0 =  8.8541878188e-12; %F/m
        A = A/4/pi/eps0./Array.n(kk).^2;     
   end
    % toc
     
    % exp_term=exp(1i*k.*R)./(R);         
    % R_sqrd_term=(1i*k.*R-1)./(R.^2);
    % 
    % r_i_per_R(:,1)=r_i(:,1)./R;
    % r_i_per_R(:,2)=r_i(:,2)./R;
    % r_i_per_R(:,3)=r_i(:,3)./R;
    
    % A_Arr = exp_term.*(k^2.*(r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)) ...
    %                   +R_sqrd_term.*(3.*r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj))); 
    % for ii=1:3     % loop over the 3 x 3 indices of the dyadic 
    %     for jj=1:3
    %         A(ii:3:end,jj:3:end) = exp(1i*k.*R)./(4*pi*R).* ...
    %                                (k^2.*(R_jk_hat{ii}.*R_jk_hat{jj}-I3(ii,jj)) ...
    %                                  +(1i*k.*R-1)./(R.^2).*(3*R_jk_hat{ii}.*R_jk_hat{jj})-I3(ii,jj));
    %     end
    % end
    
    nnzNaN=find(isnan(A));
    A((nnzNaN))=0;
    
    %% self interaction terms, including off-diagonal terms ij (in addition to ii):
    % this might be actually quite time-consuming part, since with magnetic
    % interaction we'll need to take that '\' operation quite a many times.
    % -mjh 6.7.2015
    % speed this up when particles are identical, when the time-consuming 
    % operation   squeeze(Array.alpha(ii,:,:,kk))\eye(3); needs only be done
    % once:
    % tmp=squeeze(Array.alpha(1,:,:,kk))\eye(3);
    for ii=1:Input.N  
        Aself((3*ii-2):(3*ii),(3*ii-2):(3*ii))=squeeze(Array.alpha(ii,:,:,kk))\eye(3);
    %     Aself((3*ii-2):(3*ii),(3*ii-2):(3*ii))=tmp;
    end
    [J]=find(Aself);    % find non-zero self interaction terms and input into the A(ij)
    A(J)=Aself(J);
         
    nnzNaN=find(isnan(A));
    A((nnzNaN))=0;
    if Input.Hermitian
        upperones = logical(triu(ones(size(A))));
        A(upperones) = conj(A(upperones));
    end
end