function [Array] = Make_beta_Cinfv(Array,Input)
% Make first- and second-order hyperpolarizabilities based on polarizabilities using, 
% first just assumption that the NLO processes are instantaneous, after
% which beta_ijk becomes a constant tensor. \o/

% index notations: 
% 1 = x, direction along the surface
% 2 = y, direction along the surface
% 3 = z, direction perpendicular to the surface plane 

%% SFG, the structure of the scattering nanoparticles is assumed by taking it to belong 
%  some certain symmetry, e.g. Cinf_v:  zzz, xxz=xzx=yyz=yzy, zxx=zyy  :

beta_omega3=zeros(3,3,3);   % beta is 3x3x3 tensor
C=0.5;   %1                     % for simplicity assume superbly large nonlinearity now

% handle to control whether nnz of chi is shown or not:
if isfield(Input,'NLO_modality')
    show_nnz_of_chi=true;
else
    show_nnz_of_chi=false;
end

%% this handle can be used to calculate the contribution of each separate
%  chi(2) component, since we're summing their contributions to get the
%  nonlinear source polarization and lose that information.
if isfield(Input,'beta_struct')
    switch Input.beta_struct
        case 'xxz_yyz_zxx'
            beta_omega3(1,1,3)=C;
            beta_omega3(2,2,3)=C;
            beta_omega3(3,1,1)=C;
            
        case 'xzx_yzy_zyy'
            beta_omega3(1,3,1)=C;
            beta_omega3(2,3,2)=C;
            beta_omega3(3,2,2)=C;
        case 'zzz'
            beta_omega3(3,3,3)=C;
        case 'D3h'
            beta_omega3(2,2,2)=C;                       % yyy
            beta_omega3(2,1,1)=-beta_omega3(2,2,2);     % -yxx
            beta_omega3(1,1,2)=-beta_omega3(2,2,2);     % -xxy
            beta_omega3(1,2,1)=-beta_omega3(2,2,2);     % -xyx
        case 'C1v'
%             C=10;   %1                     % for simplicity assume superbly large nonlinearity now
            % some nanoletter of L' sand T's by robert, MK et al. Fig. 3
%             beta_omega3(2,2,2)=C;                       % yyy
%             beta_omega3(2,1,1)=0.1*C;                   % -yxx
%             beta_omega3(1,1,2)=0.05*C;                  % -xxy
%             beta_omega3(1,2,1)=beta_omega3(1,1,2);      % -xyx
            
            %% 21.8. trial 
            beta_omega3(2,2,2)=0.1*C;                       % yyy
            beta_omega3(2,1,1)=C;                       % -yxx
            beta_omega3(1,1,2)=0.001*C;                       % -xxy
            beta_omega3(1,2,1)=beta_omega3(1,1,2);      % -xyx
            
            [Input.particle] = Less_is_more_particle_specs('C1v');

        otherwise
            beta_omega3(1,1,3)=C;
            beta_omega3(1,3,1)=beta_omega3(1,1,3);
            
            beta_omega3(2,2,3)=C;
            beta_omega3(2,3,2)=beta_omega3(2,2,3);
            
            beta_omega3(3,1,1)=C;
            beta_omega3(3,2,2)=beta_omega3(3,1,1);
            
            beta_omega3(3,3,3)=C;
    end
    
    if show_nnz_of_chi
        disp(['Non-zero chi components: ' num2str(nnz(beta_omega3))])
    end
    
    Array.beta=beta_omega3;
else
            beta_omega3(1,1,3)=C;
            beta_omega3(1,3,1)=beta_omega3(1,1,3);
            
            beta_omega3(2,2,3)=C;
            beta_omega3(2,3,2)=beta_omega3(2,2,3);
            
            beta_omega3(3,1,1)=C;
            beta_omega3(3,2,2)=beta_omega3(3,1,1);
            
            beta_omega3(3,3,3)=C;
            Array.beta=beta_omega3;
end

%%%%%%%%%%%% second-order hyperpolarizability i.e. gamma  %%%%
% For simplicity, we assume here an isotropic response, 
% yyzz=zzyy=zzxx=xxzz=xxyy=yyxx   (1a)
% yzyz=zyzy=zxzx=xzxz=xyxy=yxyx   (1b)
% yzzy=zyyz=zxxz=xzzx=xyyx=yxxy   (1c)
%
% xxxx=yyyy=zzzz=xxyy+xyxy+xyyx   (2)

%% CHI(3) FOR THG:
gamma_2=zeros(3,3,3,3);            % gamma is 3x3x3x3 tensor
C=1/3;                             % only 1 independent component

gamma_2(2,2,3,3)=C;
gamma_2(1,1,3,3)=gamma_2(2,2,3,3);

% should provide 21 non-zero components
for ii=1:3   
    for jj=1:3
        gamma_2(ii,ii,jj,jj)=gamma_2(1,1,3,3);  % (1)
        gamma_2(ii,jj,ii,jj)=gamma_2(1,1,3,3);  % 
        gamma_2(ii,jj,jj,ii)=gamma_2(1,1,3,3);  % 
    end
end

% we need to overwrite these still : xxxx=yyyy=zzzz=xxyy+xyxy+xyyx
for ii=1:3   
    gamma_2(ii,ii,ii,ii)=gamma_2(1,1,2,2)+gamma_2(1,2,1,2)+gamma_2(1,2,2,1); % (2)
end

Array.gamma=gamma_2;
% Array.muna=gamma_2;

%% End of codes

% plan B will be to check if more detailed approach would make any sense:
% Boyd NLO book Chapter 1:  
%   Chi(2) \propto C * Chi(1)(\omega_3) * Chi(1)(\omega_1) * Chi(1)(\omega_2)
% i.e. for molecular polarizabilities and the hyperpolarizabilities,
% similar equation should work, 
%
%   Beta_ijk \propto C * alpha_ii(omega_3) * alpha_jj(omega_1) * alpha_kk(omega_2)
%
% if omega_3 = omega_1 - omega_2 , then  alpha_kk(-omega_2) = conj(alpha_kk(omega_2))
% eV2nm=@(eV) 1239.84193./eV; 
% nm2eV=@(wl) 1239.84193./wl; 
% 
%  Alpha_lor=@(lambda_peak,wl) A0./(2*pi*c./lambda_peak-2*pi*c./wl-1i*gamma);
%  lambda_peak_vec=[500e-9 500e-9 900e-9];
%  alpha_lor_x=Alpha_lor(lambda_peak_vec(1),Array.wl);
%  alpha_lor_y=Alpha_lor(lambda_peak_vec(2),Array.wl);
%  alpha_lor_z=Alpha_lor(lambda_peak_vec(3),Array.wl);
%    
%  alpha(1,:)=alpha_lor_x;
%  alpha(2,:)=alpha_lor_y;
%  alpha(3,:)=alpha_lor_z;
%  
%  % form frequency arrays, using which the omega3 = omega1+omega2
% omega1_array=repmat(Array.wl,size(Array.wl,2),1);
% omega2_array=repmat(Array.wl.',1,size(Array.wl,2));
% omega3_array=omega1_array+omega2_array;
% 
% 
%  C=1e40;  % for simplicity assume superbly large nonlinearity now
%  % force the loop to go over only the non-zero components
%  % 
%  for ii=1:3
%      for jj=1:3
%          for kk=1:3
%                 alphaOmega3_ii=Alpha_lor(lambda_peak_vec(ii),omega3_array);
%                  
%                 beta_omega3(ii,jj,kk,:,:)=C.*alphaOmega3_ii.*(alpha(jj,:).'*alpha(kk,:));
%          end
%      end
%  end
% 
% figure
% subplot(2,2,1)
% pcolor(real(squeeze(beta_omega3(1,1,3,:,:))))
% shading interp
% colorbar
% axis square
% subplot(2,2,2)
% pcolor(real(squeeze(beta_omega3(3,1,1,:,:))))
% shading interp
% colorbar
% axis square
% subplot(2,2,3)
% pcolor(real(squeeze(beta_omega3(3,3,3,:,:))))
% shading interp
% colorbar
% axis square
% subplot(2,2,4)
% pcolor(real(squeeze(beta_omega3(1,3,1,:,:))))
% shading interp
% colorbar
% axis square

end
