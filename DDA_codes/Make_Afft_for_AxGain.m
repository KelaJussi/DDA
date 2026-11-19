function [Afft,siz]=Make_Afft_for_AxGain(kk,Array,r_i,Input)
% Calculate convolution form of Aij for single wavelength

A_BT=cell(3,3);         % BT = Block Toeplitz 
I3=eye(3,'like',r_i);   % construct eye matrix into same memory as r_i variable
eps0 =  8.8541878188e-12; %F/m
%% Define variables to speed-up calculations of A_BT:
R=sqrt(r_i(:,1).^2 + r_i(:,2).^2 + r_i(:,3).^2);
if Input.nearest_neighbour
    indexx = (abs(r_i(:,1))-Input.p(1))<1e-10;
    indexy = (abs(r_i(:,2))-Input.p(2))<1e-10;
end
R(R<1e-14)=0;   % force zero dipole R to zero, so that isinf/isnan do find the right spot after the for loop
%gain = 1e-2;
gain = Input.gain;
k = Array.k(kk);
exp_term=exp(1i*k.*R)./(R);   
%exp_term=exp(1i*k.*R)./(4*pi*R);   
R_sqrd_term=(1i*k.*R-1)./(R.^2);
gain_term = exp(gain.*k.*R);

r_i_per_R(:,1)=r_i(:,1)./R;
r_i_per_R(:,2)=r_i(:,2)./R;
r_i_per_R(:,3)=r_i(:,3)./R;

% R=sqrt(Array.r_i_conv(:,1).^2 + Array.r_i_conv(:,2).^2 + Array.r_i_conv(:,3).^2);
% R(R<1e-14)=0;   % force zero dipole R to zero, so that isinf/isnan do find the right spot after the for loop
% r_i_per_R(:,3)=Array.r_i_conv(:,3)./R;

%% Define dimensonality
switch Input.dimensionality
    case '1Dx'  
        ind_ii=1;        ind_jj=1;                         
    case '1Dy'
        ind_ii=2;        ind_jj=2;      
    case '1Dz'
        ind_ii=3;        ind_jj=3; 
    otherwise
        ind_ii=1:3;      ind_jj=1:3; 
end

%% Loop over the right 3 x 3 indices of the dyadic to calculate A_BT
for ii=ind_ii
    for jj=ind_jj
            % 26.1.2017 was like this before.. problem with the ii,jj
            % indices? -mjh s
%           A_Arr = exp_term.*(k^2.*(r_i_per_R(:,ii).*r_i_per_R(:,ii)-I3(ii,jj)) ...
%                   +R_sqrd_term.*(3.*r_i_per_R(:,ii).*r_i_per_R(:,ii)-I3(ii,jj)));  
            %omega^2*mu0 = c^2k^2mu0= k^2/eps0
            A_Arr = gain_term.*exp_term.*(k^2.*(r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)) ...
                  +R_sqrd_term.*(3.*r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)));  
          if dda_alpha_is_SI(Input,Array)   % SI polarizabilities
              A_Arr = A_Arr/4/pi./eps0./Array.n(kk).^2;

          end
          if Input.nearest_neighbour
                A_Arr(not(and(indexx,indexy))) = 0;
          end
          % find NaN location, as that's where the main dipole is and input 1/alpha into it:
          A_Arr(isnan(A_Arr))=1e-10;  % 1/1e-30 was 1e-30 before 10.3.2017
          A_Arr(isinf(A_Arr))=1e-10;  % 1/1e-30
          
          A_BT{ii,jj}=reshape(A_Arr,[Input.Nmax(1)*2-1 Input.Nmax(2)*2-1  Input.Nmax(3)*2-1]);
    end            
end
%clear A_Arr exp_term R_sqrd_term r_i_per_R

%% Check the sizes for making the correct ffts
if Input.Nmax(3)> 1
    siz.A = size(A_BT{ii});
else
    siz.A = [size(A_BT{ii,ii}) 1];
end
siz.X = [Input.Nmax(1) Input.Nmax(2) Input.Nmax(3)];
siz.total =  2.^nextpow2(siz.A + siz.X -1);
% siz.total =  siz.A + siz.X -1;

%% Fourier and inverse fourier transforms to perform the convolution:
Afft=cell(length(ind_ii),length(ind_jj));
switch Input.dimensionality
    case {'1Dx','1Dy','1Dz'}
        Afft{ii,jj}=(fftn(A_BT{ii,jj},siz.total));
    otherwise
        for jj=ind_jj
            for ii=ind_ii
                Afft{ii,jj}=(fftn(A_BT{ii,jj},siz.total));     % (jj,z)
            end
        end
end

%% plot greens function 
if 0
%plot_A(Afft{2,2})
title("yy,fft")
figure
%plot_A(A_BT{2,2})
title("yy")
figure
plot_A(Afft{1,1})
title("xx,fft")
figure
plot_A(A_BT{1,1})
title("xx")
figure
%plot_A(Afft{3,3})
title("zz,fft")
figure
%plot_A(A_BT{3,3})
title("zz")
end
end

function [Afft,siz]=Make_Afft_for_Ax_complex(k,r_i,Input)
% Calculate convolution form of Aij for single wavelength

A_BT=cell(3,3);         % BT = Block Toeplitz 
I3=eye(3,'like',r_i);   % construct eye matrix into same memory as r_i variable
eps0 =  8.8541878188e-12; %F/m
%% Define variables to speed-up calculations of A_BT:
R=sqrt(r_i(:,1).^2 + r_i(:,2).^2 + r_i(:,3).^2);
if Input.nearest_neighbour
    indexx = (abs(r_i(:,1))-Input.p(1))<1e-10;
    indexy = (abs(r_i(:,2))-Input.p(2))<1e-10;
end
R(R<1e-14)=0;   % force zero dipole R to zero, so that isinf/isnan do find the right spot after the for loop
%gain = 1e-2;
gain = Input.gain;
exp_term=exp(1i*k.*R)./(R);   
%exp_term=exp(1i*k.*R)./(4*pi*R);   
R_sqrd_term=(1i*k.*R-1)./(R.^2);
gain_term = exp(gain.*k.*R);

r_i_per_R(:,1)=r_i(:,1)./R;
r_i_per_R(:,2)=r_i(:,2)./R;
r_i_per_R(:,3)=r_i(:,3)./R;

% R=sqrt(Array.r_i_conv(:,1).^2 + Array.r_i_conv(:,2).^2 + Array.r_i_conv(:,3).^2);
% R(R<1e-14)=0;   % force zero dipole R to zero, so that isinf/isnan do find the right spot after the for loop
% r_i_per_R(:,3)=Array.r_i_conv(:,3)./R;

%% Define dimensonality
switch Input.dimensionality
    case '1Dx'  
        ind_ii=1;        ind_jj=1;                         
    case '1Dy'
        ind_ii=2;        ind_jj=2;      
    case '1Dz'
        ind_ii=3;        ind_jj=3; 
    otherwise
        ind_ii=1:3;      ind_jj=1:3; 
end

%% Loop over the right 3 x 3 indices of the dyadic to calculate A_BT
for ii=ind_ii
    for jj=ind_jj
            % 26.1.2017 was like this before.. problem with the ii,jj
            % indices? -mjh s
%           A_Arr = exp_term.*(k^2.*(r_i_per_R(:,ii).*r_i_per_R(:,ii)-I3(ii,jj)) ...
%                   +R_sqrd_term.*(3.*r_i_per_R(:,ii).*r_i_per_R(:,ii)-I3(ii,jj)));  
            %omega^2*mu0 = c^2k^2mu0= k^2/eps0
            A_Arr = gain_term.*exp_term.*(k^2.*(r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)) ...
                  +R_sqrd_term.*(3.*r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)));  
          if dda_alpha_is_SI(Input,Array)   % SI polarizabilities (was: without n^2)
              A_Arr = A_Arr/4/pi./eps0./Array.n(kk).^2;

          end
          if Input.nearest_neighbour
                A_Arr(not(and(indexx,indexy))) = 0;
          end
          % find NaN location, as that's where the main dipole is and input 1/alpha into it:
          A_Arr(isnan(A_Arr))=1e-10;  % 1/1e-30 was 1e-30 before 10.3.2017
          A_Arr(isinf(A_Arr))=1e-10;  % 1/1e-30
          
          A_BT{ii,jj}=reshape(A_Arr,[Input.Nmax(1)*2-1 Input.Nmax(2)*2-1  Input.Nmax(3)*2-1]);
    end            
end
%clear A_Arr exp_term R_sqrd_term r_i_per_R

%% Check the sizes for making the correct ffts
if Input.Nmax(3)> 1
    siz.A = size(A_BT{ii});
else
    siz.A = [size(A_BT{ii,ii}) 1];
end
siz.X = [Input.Nmax(1) Input.Nmax(2) Input.Nmax(3)];
siz.total =  2.^nextpow2(siz.A + siz.X -1);
% siz.total =  siz.A + siz.X -1;

%% Fourier and inverse fourier transforms to perform the convolution:
Afft=cell(length(ind_ii),length(ind_jj));
switch Input.dimensionality
    case {'1Dx','1Dy','1Dz'}
        Afft{ii,jj}=(fftn(A_BT{ii,jj},siz.total));
    otherwise
        for jj=ind_jj
            for ii=ind_ii
                Afft{ii,jj}=(fftn(A_BT{ii,jj},siz.total));     % (jj,z)
            end
        end
end

%% plot greens function 
if 0
%plot_A(Afft{2,2})
title("yy,fft")
figure
%plot_A(A_BT{2,2})
title("yy")
figure
plot_A(Afft{1,1})
title("xx,fft")
figure
plot_A(A_BT{1,1})
title("xx")
figure
%plot_A(Afft{3,3})
title("zz,fft")
figure
%plot_A(A_BT{3,3})
title("zz")
end
end


function plot_A(A)
    imagesc(abs(A))
end