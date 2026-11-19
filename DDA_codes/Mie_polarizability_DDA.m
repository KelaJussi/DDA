function [Array] = Mie_polarizability_DDA(Array, Input)
if nargin < 2
     r=65e-9;    % particle radius
    %--------------------------------------------------------------------------
    la=200e-9:1e-9:1600e-9;    % wavelength
    Input.range.eV = 1239.84e-9./la;
    %n_sub=1.5;
    %eh = n_sub.^2;

    eh = epsilon_BK7(Input.range.eV); 
    n_sub = sqrt(eh);
    Array.n = n_sub;

    Input.N = 1;
    Input.particle.anisotropic =0;
    Input.array_type = "";
    Input.sample.type ="Au_Drude_Mie";
    Input.dimensionality = "3D";
    plott = 1;
    plotteV  = 0; %eV or wl (nm)
elseif nargin == 2
    la = Array.wl;
    r = Input.particle.r;
    
    if Input.substrate == "BK7"
        
       eh = epsilon_BK7(Input.range.eV); 
       n_sub = sqrt(eh);
       Array.n = n_sub;
    else
        eh=Input.n^2;                      % host epsilon
        n_sub = Input.n;
        Array.n = n_sub*ones(size(Input.range.eV));
    end
    if 0
        la(la>1544.6e-9 & la<1550.8e-9) = 1544.6e-9;
    end
    c = 2.99792458e8;
   % Array.n = n_sub;
    Array.k = sqrt(eV2rads(Input.range.eV).^2/c^2.*eh);
    plott = 0;
end
    

    muh=1;                      % host mu
    mus=1;                      % particle mu
    k0=2.*pi./la;   % vacuum wave number
    %clc;
    %clear all;
    %close all
    f=0.3;     % filling factor of particles in metamaterial
    Array.alpha=zeros(Input.N,3,3,size(la,2));
    if Input.particle.anisotropic
        alphax = calc_alpha(Input.particle.rx,la);
        alphay = calc_alpha(Input.particle.ry,la);
        alphaz = calc_alpha(Input.particle.rz,la);
        for kk=1:Input.N % Input.N
            Array.alpha(kk,:,:,:) = 1e-30;
            %Array.alpha(kk,:,:,:) = 0.1*single(reshape(alphaz,1,1,1,[])).*ones(1,3,3);
            Array.alpha(kk,1,1,:) = single(alphax);
            Array.alpha(kk,2,2,:) = single(alphay);
            Array.alpha(kk,3,3,:) = single(alphaz);
            %Array.alpha(kk,1,2,:) = single(alphax);
            %Array.alpha(kk,2,1,:) = single(alphax);
            
        end
    elseif Input.dimensionality == "cyl_TM"
        Array.alpha = zeros(Input.N,1,1,size(la,2)); %[parts x, parts y, pol, wl]
        alpha = calc_alpha_cylinder(r,plott);
            for kk=1:Input.N % Input.N
                Array.alpha(kk,:,:,:) = 1e-30;
                Array.alpha(kk,1,1,:) = single(alpha);
              
            end
    elseif Input.dimensionality == "cyl_TE"
        % line dipole in the plane perpendicular to the cylinders: (p_x, p_z),
        % isotropic in that plane
        Array.alpha = zeros(Input.N,2,2,size(la,2));
        alpha = calc_alpha_cylinder(r,plott);
            for kk=1:Input.N
                Array.alpha(kk,1,1,:) = single(alpha);
                Array.alpha(kk,2,2,:) = single(alpha);
            end

    else
            if Input.array_type == "Outerwall"
               
                alpha = calc_alpha(Input,la);
                alpha_walls = calc_alpha(Input+Input.wall_increase_in_radius,la);
                 minx = min(Array.r_i(:,1));miny = min(Array.r_i(:,2));
                 maxx = max(Array.r_i(:,1));maxy = max(Array.r_i(:,2));
                 if Input.Outerwallxy ==0
                  wall_indeces = Array.r_i(:,1)==minx | Array.r_i(:,2)==miny...
                         | Array.r_i(:,1)==maxx | Array.r_i(:,2)==maxy;
                 elseif Input.Outerwallxy == 1
                    wall_indeces = Array.r_i(:,1)==minx | Array.r_i(:,1)==maxx ;
                 elseif Input.Outerwallxy == 2
                    wall_indeces = Array.r_i(:,2)==miny |  Array.r_i(:,2)==maxy;
                 end

                  x = Array.r_i(:,1);
                    y = Array.r_i(:,2);
                    B = Input.wall_thickness;
                    % Unique sorted values
                    xu = unique(x,'sorted');
                    yu = unique(y,'sorted');
                    
                    % Handle case B > number of unique values
                    Bx = min(B, numel(xu));
                    By = min(B, numel(yu));
                    
                    % B unique smallest and largest values
                    x_extreme = [xu(1:Bx); xu(end-Bx+1:end)];
                    y_extreme = [yu(1:By); yu(end-By+1:end)];
                    
                    % Logical mask
                    flag = ismember(x, x_extreme) | ismember(y, y_extreme);
                    
                                       
                    N = size(x,1);
                    
                    % --- X coordinate layers ---
                    xu = unique(x,'sorted');
                    dx = min( (1:numel(xu)).', (numel(xu):-1:1).' );
                    xLayerMap = containers.Map(xu, dx);
                    xLayer = arrayfun(@(v) xLayerMap(v), x);
                    
                    % --- Y coordinate layers ---
                    yu = unique(y,'sorted');
                    dy = min( (1:numel(yu)).', (numel(yu):-1:1).' );
                    yLayerMap = containers.Map(yu, dy);
                    yLayer = arrayfun(@(v) yLayerMap(v), y);
                    
                    % --- Final rank: proximity to any edge ---
                    rank = min(xLayer, yLayer);
                    % maxrank = max(rank);
                    % for kk=1:maxrank
                    %     rtemp = r+(maxrank-kk)/(maxrank-1)*Input.wall_increase_in_radius;
                    %     alphatemp(:,kk) = calc_alpha(Inputtemp);
                    % end
                for kk=1:Input.N % Input.N
                    Array.alpha(kk,:,:,:) = 1e-30;
                    %  Array.alpha(kk,1,1,:) = single(alphatemp(:,rank(kk)));
                    % Array.alpha(kk,2,2,:) = single(alphatemp(:,rank(kk)));
                    % Array.alpha(kk,3,3,:) = single(alphatemp(:,rank(kk)));
                     if flag(kk)
                        Array.alpha(kk,1,1,:) = single(alpha_walls);
                        Array.alpha(kk,2,2,:) = single(alpha_walls);
                        Array.alpha(kk,3,3,:) = single(alpha_walls);
                     else
                         Array.alpha(kk,1,1,:) = single(alpha);
                            Array.alpha(kk,2,2,:) = single(alpha);
                            Array.alpha(kk,3,3,:) = single(alpha);
                    end
                    
                end
            elseif Input.array_type == "Inverted_Gaussian"
                    x = Array.r_i(:,1);
                    y = Array.r_i(:,2);

                    % --- Parameters ---
                    g_min = r ;     % value at edges
                    g_max = r + Input.wall_increase_in_radius;     % value at center
                    sigma = Input.Gaussian_sigma;    % Gaussian width (adjustable)
                    

                    % Center along y
                    yc = mean(y);
                    
                    % 1D Gaussian
                    G = exp( -(y - yc).^2 / (2*sigma^2) );
                    
                    % Normalize
                    G = G / max(G);
                    
                    % Scale to desired range                    
                    Ginv = 1 - G;   % invert
                    
                    g = g_min + (g_max - g_min)*Ginv;

                    for kk = 1:Input.N
                        Array.alpha(kk,:,:,:) = 1e-30;
                        alphatemp = calc_alpha(g(kk),la);
                        Array.alpha(kk,1,1,:) = single(alphatemp);
                        Array.alpha(kk,2,2,:) = single(alphatemp);
                        Array.alpha(kk,3,3,:) = single(alphatemp);
                    end
                   abs(max(Array.alpha(:,1,1,1)));
                    abs(min(Array.alpha(:,1,1,1)));
            elseif Input.array_type == "Gaussian"
                    x = Array.r_i(:,1);
                    y = Array.r_i(:,2);

                    % --- Parameters ---
                    g_min = r ;     % value at edges
                    g_max = r + Input.wall_increase_in_radius;     % value at center
                    sigma = Input.Gaussian_sigma;    % Gaussian width (adjustable)
                    

                    % Center along y
                    yc = mean(y);
                    
                    % 1D Gaussian
                    G = exp( -(y - yc).^2 / (2*sigma^2) );
                    
                    % Normalize
                    G = G / max(G);
                    
                    % Scale to desired range                    
                   % Ginv = 1 - G;   % invert
                    
                    g = g_min + (g_max - g_min)*G;

                    for kk = 1:Input.N
                        Array.alpha(kk,:,:,:) = 1e-30;
                        alphatemp = calc_alpha(g(kk),la);
                        Array.alpha(kk,1,1,:) = single(alphatemp);
                        Array.alpha(kk,2,2,:) = single(alphatemp);
                        Array.alpha(kk,3,3,:) = single(alphatemp);
                    end
                   abs(max(Array.alpha(:,1,1,1)));
                    abs(min(Array.alpha(:,1,1,1)));
            elseif Input.array_type == "HG_wall"
                alpha = calc_alpha(Input,la);
                alpha_walls = calc_alpha(Input+Input.wall_increase_in_radius,la);
                xmode=Input.HG_wall_mode(1);ymode=Input.HG_wall_mode(2);
                % maxx = max(Array.r_i(:,1));maxy = max(Array.r_i(:,2));
                if xmode
                    xindeces = ceil((Input.arrayN(1)+1).*(1:xmode)./(xmode+1));
                else
                    xindeces = NaN;
                end
                if ymode
                  yindeces = ceil((Input.arrayN(2)+1).*(1:ymode)./(ymode+1));
                else
                    yindeces = NaN;
                end
                

                for kk=1:Input.N % Input.N
                    Array.alpha(kk,:,:,:) = 1e-30;
                     if any(mod(kk,Input.arrayN(1)) == xindeces)  || any(ceil(kk/Input.arrayN(1))==yindeces)
                        Array.alpha(kk,1,1,:) = single(alpha_walls);
                        Array.alpha(kk,2,2,:) = single(alpha_walls);
                        Array.alpha(kk,3,3,:) = single(alpha_walls);
                     else
                         Array.alpha(kk,1,1,:) = single(alpha);
                            Array.alpha(kk,2,2,:) = single(alpha);
                            Array.alpha(kk,3,3,:) = single(alpha);
                    end
                    
                end


            else
                alpha = calc_alpha(Input,la,plott);
                for kk=1:Input.N % Input.N
                    Array.alpha(kk,:,:,:) = 1e-30;
                    Array.alpha(kk,1,1,:) = single(alpha);
                    Array.alpha(kk,2,2,:) = single(alpha);
                    Array.alpha(kk,3,3,:) = single(alpha);
                end
                
            end
    end
    %alphajtk = alphajtk.*V;
    
    %% Alpha
   
    
    
    eps0 =  8.8541878188e-12;

    n_sub_r = reshape(n_sub, 1, 1, 1, []);   % broadcast over the wavelength dimension
    Array.alpha = eps0*n_sub_r.*Array.alpha;
    if 0 && not(ishandle(303))
        figure(303)
        plot(la,real(squeeze(alpha*eps0)),la,imag(eps0*alpha))
    end
   
    


    function alpha = calc_alpha_cylinder(a, plott)
        if  Input.sample.type=="Au_Mie"
            % Johnson & Christy gold data, next to this file in DDA_codes
            filename = fullfile(fileparts(mfilename('fullpath')), 'Johnson_gold.csv');
            if ~isfile(filename)
                error('Johnson_gold.csv not found next to %s.', mfilename);
            end
            data = load(filename);
            lambda = data(:,1)*1e-6;
            n = data(:,2); kappa = data(:,4);
            expE1 = (n.^2-kappa.^2);
            expE2 = (2*n.*kappa);
            exp_eps = expE1 + 1i*expE2;                         %complex permittivity data points
            %omega_vec = h*c./lambda;
            
            exp_eps_vec = interp1(lambda,exp_eps,la);
            es = exp_eps_vec;
        elseif Input.sample.type == "glass_Mie"
            kappa = [Input.particle_gain,Input.particle_gain];
            n = [Input.Mie_glass_n,Input.Mie_glass_n];
            lambda =[0,1e-5];
            expE1 = (n.^2-kappa.^2);
            expE2 = (2*n.*kappa);
            exp_eps = expE1 + 1i*expE2;                         %complex permittivity data points
            %omega_vec = h*c./lambda;
            
            exp_eps_vec = interp1(lambda,exp_eps,real(la));
            es = exp_eps_vec;

        elseif Input.sample.type == "Au_Drude_Mie"
            omega = Input.range.eV;
            %omega = 1240e-9/la+Input.gain; % eV
            omega_p = 9; %eV
            gamma = 70e-3; %eV, Berciaud, http://www.wave-scattering.com/drudefit.html
            eps_inf = 9; %common according to copilot
             es = eps_inf - omega_p.^2./(omega.^2+1i*gamma*omega);
             n = sqrt(abs(es)+real(es)/2);
             kappa = sqrt(abs(es)-real(es)/2);
        end
            kh = sqrt(eh).*k0;
            x = sqrt(eh) .* k0 * a;
            m = sqrt(es ./ eh);

            % Mie theory of an infinite cylinder, lowest term of each polarization.
            % Scattered field = s_n * H_n(kr) e^{i n phi} per unit incident amplitude:
            %   cyl_TM (E along the axis), n = 0:
            %     s0 = [J0'(x)J0(mx) - m J0(x)J0'(mx)] / [m H0(x)J0'(mx) - H0'(x)J0(mx)]
            %     alpha = -4i*s0*eps0*n^2/k^2      -> pi a^2 eps0 (es - eh) for small a
            %   cyl_TE (E perpendicular to the axis), n = +-1:
            %     s1 = [m J1'(x)J1(mx) - J1(x)J1'(mx)] / [H1(x)J1'(mx) - m H1'(x)J1(mx)]
            %     alpha = -8i*s1*eps0*n^2/k^2      -> 2 pi a^2 eps0 eh (es - eh)/(es + eh)
            % Both include the radiation reaction, as the matrix builders expect.
            % The values returned here are alpha/(eps0*n): the factor eps0*n is
            % applied to every Mie polarizability at the end of this file.
            legacy = isfield(Input.particle,'cyl_alpha_legacy') && Input.particle.cyl_alpha_legacy;
            if ~legacy || Input.dimensionality == "cyl_TE"
                dJ = @(n,z) (besselj(n-1,z) - besselj(n+1,z))/2;
                dH = @(n,z) (besselh(n-1,1,z) - besselh(n+1,1,z))/2;
                if Input.dimensionality == "cyl_TE"
                    s1 = (m.*dJ(1,x).*besselj(1,m.*x) - besselj(1,x).*dJ(1,m.*x)) ./ ...
                         (besselh(1,1,x).*dJ(1,m.*x) - m.*dH(1,x).*besselj(1,m.*x));
                    alpha = -8i*s1.*sqrt(eh)./kh.^2;
                else
                    s0 = (dJ(0,x).*besselj(0,m.*x) - m.*besselj(0,x).*dJ(0,m.*x)) ./ ...
                         (m.*besselh(0,1,x).*dJ(0,m.*x) - dH(0,x).*besselj(0,m.*x));
                    alpha = -4i*s0.*sqrt(eh)./kh.^2;
                end
                if plott
                    plot(2*pi./k0*1e9,real(alpha),2*pi./k0*1e9,imag(alpha))
                end
                return
            end
            % Input.particle.cyl_alpha_legacy = true: the cyl_TM expression used
            % before. It is too small by the factor (n_cyl + n_host) for thin
            % cylinders (alpha -> pi a^2 eps0 n_h (n_s - n_h)); kept only to
            % reproduce old results.
            
            J0x   = besselj(0,x);
            J0mx  = besselj(0,m.*x);
            
            J0x_p  = -besselj(1,x);        % derivative
            J0mx_p = -besselj(1,m.*x).*m;
            
            H0x   = besselh(0,1,x);
            H0x_p = -besselh(1,1,x);
            
            num = m .* J0mx.*J0x_p - J0x .* J0mx_p;
            den = m .* J0mx .* H0x_p - H0x .* J0mx_p;
            
            alpha = (4i/1) * (num ./ den)./kh.^2;

            if plott

                plot(2*pi./k0*1e9,real(alpha),2*pi./k0*1e9,imag(alpha))
            end
    end

end