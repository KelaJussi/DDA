function [Array] = MLWA_polarizability_DDA(Array, Input)
if nargin < 2
     r=30*10^-9;    % particle radius
    %--------------------------------------------------------------------------
    la=200e-9:1e-9:1200e-9;    % wavelength

    n_sub=1.5;
    eh = n_sub.^2;
    Input.N = 1;
elseif nargin == 2
    la = Array.wl;
    r = Input.particle.r;
    
    if Input.substrate == "BK7"
        
       eh = epsilon_BK7(Input.range.eV); 
       n_sub = sqrt(eh);
    else
        eh=Input.n^2;                      % host epsilon
        n_sub = Input.n;
    end
    if 0
        la(la>1544.6e-9 & la<1550.8e-9) = 1544.6e-9;
    end
end
    c = 2.99792458e8;
    Array.k = sqrt(eV2rads(Input.range.eV).^2/c^2.*eh);
    muh=1;                      % host mu
    mus=1;                      % particle mu
    k0=2.*pi./la;   % vacuum wave number
    %clc;
    %clear all;
    %close all
    f=0.3;     % filling factor of particles in metamaterial
    Array.alpha=zeros(Input.N,3,3,size(la,2),'single');
    if Input.particle.anisotropic
        alphax = calc_alpha(Input.particle.rx);
        alphay = calc_alpha(Input.particle.ry);
        alphaz = calc_alpha(Input.particle.rz);
        for kk=1:Input.N % Input.N
            Array.alpha(kk,:,:,:) = 1e-30;
            %Array.alpha(kk,:,:,:) = 0.1*single(reshape(alphaz,1,1,1,[])).*ones(1,3,3);
            Array.alpha(kk,1,1,:) = single(alphax);
            Array.alpha(kk,2,2,:) = single(alphay);
            Array.alpha(kk,3,3,:) = single(alphaz);
            %Array.alpha(kk,1,2,:) = single(alphax);
            %Array.alpha(kk,2,1,:) = single(alphax);
            
        end
    else
        alpha = calc_alpha(r);
        for kk=1:Input.N % Input.N
            Array.alpha(kk,:,:,:) = 1e-30;
            Array.alpha(kk,1,1,:) = single(alpha);
            Array.alpha(kk,2,2,:) = single(alpha);
            Array.alpha(kk,3,3,:) = single(alpha);
        end
    end
    %alphajtk = alphajtk.*V;
    
    %% Alpha
   
    
    
    eps0 =  8.8541878188e-12;
    Array.alpha = eps0*Array.alpha;
    if 0 && not(ishandle(303))
        figure(303)
        plot(la,real(squeeze(alpha*eps0)),la,imag(eps0*alpha))
    end
   
    function [alpha] = calc_alpha(r)

        if Input.sample.type == "Au_Drude_MLWA"
            omega = Input.range.eV;
            %omega = 1240e-9/la+Input.gain; % eV
            omega_p = 9; %eV
            gamma = 70e-3; %eV, Berciaud, http://www.wave-scattering.com/drudefit.html
            eps_inf = 9; %common according to copilot
             es = eps_inf - omega_p.^2./(omega.^2+1i*gamma*omega);
             n = sqrt(abs(es)+real(es)/2);
             kappa = sqrt(abs(es)-real(es)/2);
        end
        alpha = 4*pi*r/3*(es-eh)./(eh+(es-eh)/3);

        
        
    
    end
end