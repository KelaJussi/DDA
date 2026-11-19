function [alpha] = calc_alpha(Input, la,plott)
        if nargin ==2
            plott = 0;
        end
        if Input.substrate == "BK7"
        
           eh = epsilon_BK7(Input.range.eV); 
           n_sub = sqrt(eh);
           Array.n = n_sub;
        else
            eh=Input.n^2;                      % host epsilon
            n_sub = Input.n;
            Array.n = n_sub*ones(size(Input.range.eV));
        end
        muh = 1;mus = 1;
        f=0.3;
        r = Input.particle.r;
        k0=2.*pi./la;   % vacuum wave number
        
        W=2*pi*3*10^8./la;          % angular freq.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        e=1.6e-19;
        h=6.62e-34;
    
          %________Au__________________________________
        
        Wp=2*pi*9.03*e/h;
        Vf=1.39e6;
        
        f0=0.760;
        f1=0.024;
        f2=0.010;
        f3=0.071;
        f4=0.601;
        f5=4.384;
        
        G0=2*pi*0.053*e/h+(Vf/((4/3)*r));
        G1=2*pi*0.241*e/h;
        G2=2*pi*0.345*e/h;
        G3=2*pi*0.870*e/h;
        G4=2*pi*2.494*e/h;
        G5=2*pi*2.214*e/h;
        
        W1=2*pi*0.415*e/h;
        W2=2*pi*0.830*e/h;
        W3=2*pi*2.969*e/h;
        W4=2*pi*4.304*e/h;
        W5=2*pi*13.32*e/h;
        
       esf=1-(     (f0*(Wp^2))./(W.*(W+i*G0))     );
    
        E1=f1*(Wp^2)./(W1^2-W.^2-i*W*G1);
        E2=f2*(Wp^2)./(W2^2-W.^2-i*W*G2);
        E3=f3*(Wp^2)./(W3^2-W.^2-i*W*G3);
        E4=f4*(Wp^2)./(W4^2-W.^2-i*W*G4);
        E5=f5*(Wp^2)./(W5^2-W.^2-i*W*G5);
        
        esb=E1+E2+E3+E4+E5;
        
        es= esf+esb; %Epsilon, sphere
        
        %% exp data from Johnson
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
            
            exp_eps_vec = interp1(lambda,exp_eps,la);
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


        
        
    
        
        x=sqrt(eh.*muh).*k0.*r;
        %x=sqrt(muh).*k0.*r;
        xs=sqrt(es.*mus).*k0.*r;
        %--------------------------------------------------------------------------
        
        aa=sqrt(pi./(2*xs)).*besselj(1.5,xs).*sqrt(pi./(2*x)).*(x.*besselj(0.5,x)-besselj(1.5,x));
        bb=sqrt(pi./(2*x)).*besselj(1.5,x).*sqrt(pi./(2*xs)).*(xs.*besselj(0.5,xs)-besselj(1.5,xs));
        cc=sqrt(pi./(2*xs)).*besselj(1.5,xs).*sqrt(pi./(2*x)).*(x.*besselh(0.5,1,x)-besselh(1.5,1,x));
        dd=sqrt(pi./(2*x)).*besselh(1.5,1,x).*sqrt(pi./(2*xs)).*(xs.*besselj(0.5,xs)-besselj(1.5,xs));
        %% test JTK, Failed
        m = sqrt(es.*mus./(eh.*muh));
        j1x=sqrt(pi./(2.*x)).*besselj(1.5,x); j1xs=sqrt(pi./(2.*xs)).*besselj(1.5,xs);
        j1mx=sqrt(pi./(2.*m.*x)).*besselj(1.5,m.*x);
        y1x = sqrt(pi./(2.*x)).*bessely(1.5,x);
        h1x = j1x+1i.*y1x;
    
        a1 = (muh.*m.^2.*j1mx.*conj(x.*j1x)-mus.*j1x.*conj(m.*x.*j1mx))./...
            (mus.*m.^2.*j1mx.*conj(x.*h1x)-muh.*h1x.*conj(m.*x.*j1mx));
    
        %alphajtk = 1i*3/2*r.^3./x.^3.*a1;                       %failed test
        
        %% working
    %     nmax = 2;
    %     for n=1:nmax; 
    %                 aa=aa + sqrt(pi./(2.*xs.*x)).*besselj(n+0.5,xs).*(x.*besselj(n-0.5,x)-n.*besselj(n+0.5,x));
    %                 bb=bb + sqrt(pi./(2.*xs.*x)).*besselj(n+0.5,x).*(xs.*besselj(n-0.5,xs)-n.*besselj(n+0.5,xs));
    %                 cc=cc + sqrt(pi./(2.*xs.*x)).*besselj(n+0.5,xs).*(x.*besselh(n-0.5,1,x)-n.*besselh(n+0.5,1,x));
    %                 dd=dd + sqrt(pi./(2.*xs.*x)).*besselh(n+0.5,1,x).*(xs.*besselj(n-0.5,xs)-n.*besselj(n+0.5,xs));
    %     end
    
        TE=(aa.*es-bb.*eh)./(cc.*es-dd.*eh);
        TH=(aa.*mus-bb.*muh)./(cc.*mus-dd.*muh);
        
        alfa=(3/2)*1i*(TE./(x.^3)); % electric polarizability
        %alfa=(9/2)*1i*(TE./(x.^3)); % electric polarizability
    
        V = 4*pi*r^3.*n_sub;
        alpha = alfa.*V;

        %alpha = (3/2)*1i * TE ./ (k0.^3 .* (eh.*muh).^(3/2));
        

    alphaa = reshape(alpha, 1,1,[]);
    
    
    zero = zeros(1,1,length(alpha));
    
    %alpha = [alphaa zero zero;zero alphaa zero;zero zero alphaa];
    
    %Part(rr).alpha = alphaa;
    alphaH=(3/2)*i*(TH./(x.^3));  % magnetic polarizability
    %alfa=(es-eh)./(es+2*eh);
    
    Qsca=(6./(x.^2)).*(abs(TE).^2+abs(TH).^2);   % scattering coefficient 
    Qext=(6./(x.^2)).*real(TE);                  % extinction coefficient 
    Qext2=(6./(x.^2)).*real(TE+TH);                  % extinction coefficient
    Qabs=Qext-Qsca;                                     % absorption coefficient 
    
    eff=eh.*(x.^3+3*i*f.*TE)./(x.^3-1.5*i*f.*TE);       % effective electric polarizability
    mueff=muh.*(x.^3+3*i*f.*TH)./(x.^3-1.5*i*f.*TH);     % effective magnetic polarizability
    reff=real(eff);
    ieff=imag(eff);
    rmueff=real(mueff);
    imueff=imag(mueff);
    
    eta=sqrt(mueff./eff);
    realeta=real(eta);
    imageta=imag(eta);
    abseta=abs(eta);
    
    ref=(eta-1)./(eta+1);
    absref=abs(ref);
    
    abse=abs(eff);
    ange=angle(eff);
    absh=abs(mueff);
    angh=angle(mueff);
    nt=sqrt(abse.*absh).*exp(i.*(ange+angh)/2);
    
    
    %plot(la,real(es),la,imag(es),'r--');
    %plot(la,real(nn),la,imag(nn),'r--');
    %plot(la,reff,la,ieff,'r--');
    if plott ==1
        figure
        %plot(la,real(alpha),la,imag(alpha));
        ext = k0.*reshape(imag(alphaa),1,[]);
        if plotteV
            plot(Input.range.eV,ext./(pi.*r.^2))
            hold on
            plot(Input.range.eV,Qext,Input.range.eV,Qext2)
            xlabel("freq (eV)")
        else
            plot(la*1e9,ext./(pi.*r.^2))
            hold on
            plot(la*1e9,Qext,la*1e9,Qext2)
            xlabel("Wl (nm)")
        end
        legend(["From alpha","Qext with magnetic", "Qext, no magnetic"])
    end
    %hold on
     %plot(la,real(alphajtk),la,imag(alphajtk));
     %legend(["old real", "old, imag", "new real", "new imag"])
    %plot(la,reff,la,ieff,'--');
    %plot(la,real(alpha),la,imag(alpha));
     %plot(la,real(es),la,imag(es),'r--');
    %plot(la,Qext,la,Qabs,'r--');
end
