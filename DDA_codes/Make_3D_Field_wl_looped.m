function [Field]=Make_3D_Field_wl_looped(Input,Array,jj)
% make 3Nx1 incident field array for given wavelength (jj)

Save_function_workspace_when_debugging()

% transform input degs into radians
theta_Rad=degtorad(Input.theta);
phi_Rad=degtorad(Input.phi);
psi_Rad=degtorad(Input.psi);
Field.wl = Array.wl(jj);
%% Make plane wave excitation:
Field.k = [ cos(phi_Rad).*sin(theta_Rad) sin(phi_Rad).*sin(theta_Rad) -cos(theta_Rad)];       % 3D input wave's unit k-vector

% Form 3D rotation matrix along the k-vector: [w.g. wikipedia: Rotation matrix]
psi = psi_Rad;  
R_rot = [cos(psi)+Field.k(1)^2*(1-cos(psi))...
         Field.k(1)*Field.k(2)*(1-cos(psi))-Field.k(3)*sin(psi)...
         Field.k(1)*Field.k(3)*(1-cos(psi))+Field.k(2)*sin(psi);...
         Field.k(2)*Field.k(1)*(1-cos(psi))+Field.k(3)*sin(psi)...
         cos(psi)+Field.k(2)^2*(1-cos(psi))...
         Field.k(3)*Field.k(1)*(1-cos(psi))-Field.k(2)*sin(psi);...
         Field.k(3)*Field.k(2)*(1-cos(psi))-Field.k(1)*sin(psi)...
         Field.k(3)*Field.k(2)*(1-cos(psi))+Field.k(1)*sin(psi)...
         cos(psi)+Field.k(3)^2*(1-cos(psi))];

% Make incident field transverse to k, by taking cross product of k and z-axis:
z_axis=[0 0 1];                     % z-axis
E0_trans=cross(Field.k,z_axis)';    % Make E0 transverse to k and z-axis
 
if nnz(E0_trans) == 0
    % make incident plane, which depends on phi
    phi = phi_Rad;
    Inc_plane=[ cos(phi) -sin(phi) ; sin(phi) cos(phi) ]*[1 ; 0];
    E0_trans=cross(Field.k,[Inc_plane ; 0])';    % Make E0 transverse to k and y-axis
end
E0_trans=E0_trans./norm(E0_trans);  % normalize E0
Field.pol = R_rot*E0_trans;         % Final Input field's polarization
switch Input.polarization
    case "RC"
        Field.pol(2) = Field.pol(2)*1i;
    case "linear"
        Field.pol(2) = Field.pol(2)*1;
    case "LC"
        Field.pol(2) = -Field.pol(2)*1i;
end
%% build E-field into correct N x 3 x wl_n format:
% let's forget the -iwt relation, and play with ~time-harmonic excitation

% k.r
Array_in = Array;                 % with the z positions (layered background)
Array.r_i(:,3)=0;
K=repmat(Field.k.*Array.k(jj)./Array.n(jj),size(Array.r_i,1),1);
K_dot_R = sum(K.*Array.r_i,2);

Field_Phase=exp(1i.*K_dot_R);            

% %% build E-field into correct 3N x 1 x wl_n format:
% Field.Einc_3N=zeros(3*Input.N,1,1,'single');
% 
% for pol=1:3
%     ind = pol:3:(3*Input.N);
%     Field.Einc_3N(ind,1,:)=single(Field.pol(pol).*Field_Phase);
% end


%% build E-field into correct 3N x 1 x wl_n (or N x wl_n) format:
Field.Einc=zeros(Input.N,1,1,'like',Array.r_i);

switch Input.dimensionality
    case '1Dx'
        Field.Einc(:,1,1)=single(Field.pol(1).*Field_Phase);     %  Ex
    case '1Dy'
        Field.Einc(:,1,1)=single(Field.pol(2).*Field_Phase);     %  Ey
    case '1Dz'
        Field.Einc(:,1,1)=single(Field.pol(3).*Field_Phase);     %  Ez
    case 'cyl_TM'
        Field.Einc(:,1,1)=single(Field.pol(2).*Field_Phase);     %  Ey
        L = dda_layers(Input, Array_in);
        if L.enable
            % plane wave from the top half-space: field inside the host layer
            % (transmission and all reflections of the stack, s polarization)
            F = dda_layer_field(Array.k(jj)/Array.n(jj), Array.n(jj), L, ...
                                Field.k(1)*Array.k(jj)/Array.n(jj), L.z, 'TM');
            Field.Einc(:,1,1)=single(Field.pol(2).*Field_Phase.*F.y);
        end
    case 'cyl_TE'
        % line dipoles in the xz plane, interleaved [Ex1 Ez1 Ex2 Ez2 ...]
        Field.Einc=zeros(2*Input.N,1,1,'like',Array.r_i);
        Field.Einc(1:2:end)=single(Field.pol(1).*Field_Phase);   %  Ex
        Field.Einc(2:2:end)=single(Field.pol(3).*Field_Phase);   %  Ez
        L = dda_layers(Input, Array_in);
        if L.enable
            % p-polarized part of the plane wave inside the host layer
            F = dda_layer_field(Array.k(jj)/Array.n(jj), Array.n(jj), L, ...
                                Field.k(1)*Array.k(jj)/Array.n(jj), L.z, 'TE');
            a_p = -(Field.pol(1)*(-Field.k(3)) + Field.pol(3)*Field.k(1));
            Field.Einc(1:2:end)=single(a_p.*Field_Phase.*F.x);
            Field.Einc(2:2:end)=single(a_p.*Field_Phase.*F.z);
        end
    case '2D'
        Field.Einc=zeros(2*Input.N,1,1,'like',Array.r_i);
        for pol=1:2
            Field.Einc(pol:2:end)=single(Field.pol(pol).*Field_Phase); %if reshaped: pol,x,y, is correct for dimensions
        end
    case '3D'
        Field.Einc=zeros(3*Input.N,1,1,'like',Array.r_i);
        for pol=1:3
            Field.Einc(pol:3:end)=single(Field.pol(pol).*Field_Phase); %if reshaped: pol,x,y, is correct for dimensions
        end
    otherwise
        Field.Einc=zeros(3*Input.N,1,1,'like',Array.r_i);
        for pol=1:3
            Field.Einc(pol:3:end)=single(Field.pol(pol).*Field_Phase);
        end
end

%% normalize incident field component:
switch Input.dimensionality
    case {'1Dx' , '1Dy', '1Dz'}
        if jj <2
            disp('Excitation field component normalized to the excitation component.')
        end
        %Field.Einc=Field.Einc./max(Field.Einc);
        Field.normalization_scaler=max(Field.Einc);
%         Field.Einc=Field.Einc./max(abs(Field.Einc));
%         Field.normalization_scaler=max(abs(Field.Einc));
end

%% Make predifined structured incident field
xsize = Input.arrayN(1);ysize = Input.arrayN(2);
if Input.dimensionality=="3D"
    dim = 3;
elseif Input.dimensionality=="2D" || Input.dimensionality=="cyl_TE"
    dim = 2;
else
    dim = 1;
end
switch Input.preset_field
    case 'none'
        ;
    case 'PC_mirror_reflection_test'
        Field.Einc = reshape(Field.Einc,dim,xsize,ysize);
        n_1 = (ysize-Input.PC_mirror_thickness) / 2 + 1;
        if n_1<0
            n_1 = 2;
        end
        Field.Einc(:,:,n_1:end) = 0;
        Field.Einc = reshape(Field.Einc,[],1);

    case'sine_y'
        order = Input.preset_sine_order;
        Field.Einc = reshape(Field.Einc,dim,xsize,ysize);
        xvec = [1:xsize]/xsize-0.5;
        Field.Einc = Field.Einc.*reshape(cos(xvec*order*pi+Input.preset_sine_phase),1,[],1);
        if 0
            figure(101)
            imagesc(reshape(-Field.Einc(1,:,:),xsize,ysize))
            colorbar
        end
        Field.Einc = reshape(Field.Einc,[],1);
    case'exp_y'
        order = Input.preset_sine_order;
        Field.Einc = reshape(Field.Einc,dim,xsize,ysize);
        Field.Einc = Field.Einc.*reshape(exp(1i*[1:xsize]/xsize*order*pi),1,[],1);
        if 0
            figure(101)
            imagesc(reshape(abs(-Field.Einc(1,:,:)),xsize,ysize))
            colorbar
        end
        Field.Einc = reshape(Field.Einc,[],1);
    case 'inverse_gaussian_y'
        vec = [1:xsize]/xsize; %normalized, 0 to 1
        Field.Einc = reshape(Field.Einc,dim,xsize,ysize);
        sigma2 = 0.05; mu = 0.5; 
        Field.Einc = Field.Einc.*reshape(0.5-exp(-1*(vec-mu).^2./2/sigma2),1,[],1); %./sqrt(2*pi*sigma2)
        if 0
            figure(101)
            imagesc(reshape(-Field.Einc(1,:,:),xsize,ysize))
            colorbar
        end
        Field.Einc = reshape(Field.Einc,[],1);
    case 'gaussian_y' %Gaussian only in magnitude, not in phase 
        vec = [1:xsize]/xsize; %normalized, 0 to 1
        Field.Einc = reshape(Field.Einc,dim,xsize,ysize);
        sigma2 = 1; mu = 0.5; Rz = 0.02;
        Field.Einc = Field.Einc.*reshape(exp((vec-mu).^2.*(-1./2/sigma2+1i/2/Rz)),1,[],1); %./sqrt(2*pi*sigma2)
        %Field.Einc = Field.Einc.*exp(-1i*(vec-mu).^2)
        if 0
            figure(101)
            %imagesc(reshape(-Field.Einc(1,:,:),xsize,ysize))
            Einc = reshape(Field.Einc(1,:,:),xsize,ysize);
            subplot(121)
            imagesc(abs(Einc));colorbar
            title(['Incident field, pol: ', num2str(1), ', absolute value'])
            subplot(122)
            imagesc(angle(Einc)/pi);colorbar
             title("Incident field, phase/pi")
        end
        Field.Einc = reshape(Field.Einc,[],1);
    case 'hermite_gaussian'
        % Field profile always the size of the array
        %x and y are all over the place
        w0 =0.5; %2x the standard deviation of the gaussian
        w0x= w0*Input.Hermite_scale_x;w0y =w0*Input.Hermite_scale_y;
        vecx = reshape([1:xsize]/xsize*2-1+Input.Hermite_offset_x,1,[]); %normalized, -1 to 1
        vecy =reshape([1:ysize].'/ysize*2-1+Input.Hermite_offset_y,1,1,[]); %normalized, -1 to 1
        Field.Einc = reshape(Field.Einc,dim,xsize,ysize);
       % x0=linspace(-2,2,501);
       % [X,Y]=meshgrid(x0,x0); %create a grid of width 4 units x 4 units with 501x501 mesh
        E0=exp(-(vecx.^2/w0x^2+vecy.^2/w0y^2).*(1+1i./Input.Hermite_R)); %create the gaussian
        

        Ax = hermiteH(Input.Hermite_order_x,sqrt(2)*vecx/w0x); 
        % Ax = reshape(Ax(:,:,:,1),1,1,[]);
        Ax = Ax(:,:,:,1);

        vecy = reshape(vecy,1,[]);
        Ay= hermiteH(Input.Hermite_order_y,sqrt(2)*vecy/w0y);
        Ay = reshape(Ay(:,:,1,1),1,1,[]);
       Ay = Ay(:,:,:,1);
        Field.Einc=Field.Einc.*Ay.*Ax.*E0; 
        if 0
            figure(102)
            subplot(121)
            Einc = reshape(Field.Einc(1,:,:),xsize,ysize);
            imagesc(abs(Einc).^2); colorbar;
            title(['Incident field, pol: ', num2str(1), ', absolute value'])
                subplot(122)
                imagesc(angle(Einc)/pi);colorbar
                 title("Incident field, phase/pi")
                 
        end
        Field.Einc = reshape(Field.Einc,[],1);
    case 'gaussian_proper'
        n = 1;
        w0 = Input.Gaussian_w0;
        z = Input.Gaussian_z;
        wl = 2*pi./Array.k(jj);
        zR = pi*w0^2*n/wl;
        wz = w0*sqrt(1+(z/zR).^2);
        Rz = z*(1+(zR/z).^2);
        if isnan(Rz)
            Rz = inf;
        end
        phiz = atan(z/zR); %Guoy phase; is const over r, not included
        % vecx = reshape([1:xsize]*Input.p(1)-xsize*Input.p(1)/2,1,1,[]); %true coordinates, 0 at center
        % vecy = reshape([1:ysize]*Input.p(2)-ysize*Input.p(2)/2,1,[]); %true coordinates, 0 at centerend
        % Field.Einc = reshape(Field.Einc,dim,ysize,xsize);
        % E0 = w0/wz*exp(-(vecx.^2+vecy.^2).*(1./wz.^2+1i*Array.k(jj)./(2*Rz)));
        % 
        % Ax = hermiteH(Input.Hermite_order_x,sqrt(2)*vecx/wz); 
        %  Ax = reshape(Ax(:,:,1,1),1,1,[]);
        % Ax = Ax(:,:,:,1);
        % 
        % Ay= hermiteH(Input.Hermite_order_y,sqrt(2)*vecy/wz);
        % Ay = Ay(:,:,:,1);
        % Field.Einc=Field.Einc.*Ay.*Ax.*E0; 
        % Field.Einc = Field.Einc.*E0;
        % Field.Einc = reshape(Field.Einc,[],1);


        %% second try, try different dimensions
        vecx = reshape([1:xsize]*Input.p(1)-xsize*Input.p(1)/2,1,[]); %true coordinates, 0 at center
        vecy = reshape([1:ysize]*Input.p(2)-ysize*Input.p(2)/2,1,1,[]); %true coordinates, 0 at centerend
        Field.Einc = reshape(Field.Einc,dim,xsize,ysize);
        E0 = w0/wz*exp(-(vecx.^2+vecy.^2).*(1./wz.^2+1i*Array.k(jj)./(2*Rz)));

        Ax = hermiteH(Input.Hermite_order_x,sqrt(2)*vecx/wz); 
         %
        Ax = Ax(:,:,:,1);
        
        Ay= hermiteH(Input.Hermite_order_y,sqrt(2)*reshape(vecy,1,[])/wz);
        Ay = reshape(Ay(:,:,1,1),1,1,[]);
        Ay = Ay(:,:,:,1);
        Field.Einc=Field.Einc.*Ay.*Ax.*E0; 
        %Field.Einc = Field.Einc.*E0;
        Field.Einc = reshape(Field.Einc,[],1); %first goes through pol, then x, then y
    case 'hermite_JTK_tests'
       % Field profile always the size of the array
        %x and y are all over the place
        w0 =0.5; %2x the standard deviation of the gaussian
        w0x= w0*Input.Hermite_scale_x;w0y =w0*Input.Hermite_scale_y;
        vecx = reshape([1:xsize]/xsize*2-1+Input.Hermite_offset_x,1,[]); %normalized, -1 to 1
        vecy =reshape([1:ysize].'/ysize*2-1+Input.Hermite_offset_y,1,1,[]); %normalized, -1 to 1
        Field.Einc = reshape(Field.Einc,dim,xsize,ysize);
       % x0=linspace(-2,2,501);
       % [X,Y]=meshgrid(x0,x0); %create a grid of width 4 units x 4 units with 501x501 mesh
        E0=exp(-(vecx.^2/w0x^2+vecy.^2/w0y^2).*(1+1i./Input.Hermite_R)); %create the gaussian
        
        hermite_phase = ones(size(vecy)); 
        %hermite_phase(1:ceil(length(hermite_phase)/2))=-1;
        hermite_phase(ceil(length(hermite_phase)/2-5):ceil(length(hermite_phase)/2+5))=-1;
        figure
        %imagesc(squeeze(hermite_phase))
        fourE0 = fftn(squeeze(E0.*hermite_phase),[1024,1024]);
        fourE0 = fftshift(fourE0);
        imagesc(squeeze(abs(fourE0)))
        E0modulated = ifftn(hermite_phase.*fourE0);
         Field.Einc=Field.Einc.*E0modulated; 
         if 1
            figure(101)
            subplot(121)
            Einc = reshape(Field.Einc(1,:,:),xsize,ysize);
            imagesc(abs(Einc).^2); colorbar;
            title(['Incident field, pol: ', num2str(1), ', absolute value'])
                subplot(122)
                imagesc(angle(Einc)/pi);colorbar
                 title("Incident field, phase/pi")
                 return
        end

    case 'nullfield'
        Field.Einc = 0*Field.Einc;
end
function a=hermiteH(n,x)
            %generate the nth hermite polynomial of x
            %pratheeb N.P (2025). hermite gaussian beam (https://www.mathworks.com/matlabcentral/fileexchange/71642-hermite-gaussian-beam), MATLAB Central File Exchange. Retrieved February 12, 2025. 
            % Produces a copy of the array at fourth dimension when n>1,
            % for no apparent reason.
            m=0:floor(n/2);
            [q1,q2]=size(m);
            s=ndims(x);
            [g1,g2]=size(x);
            X=repmat(x,[ones(1,s),q1, q2]);
            m=permute(repmat(m,[g1,1,g2]),[1,3,2]);
            a=factorial(n)*sum((-1).^m./(factorial(m).*factorial(n-2*m)).*(2*X).^(n-2*m),3);
        end 
%% remove unnecessary fieds from Array struct:
% Field=rmfield(Field, 'Einc_Nx2');
% Field=rmfield(Field, 'Phase');


% Function to save the function workspace when debugging
function Save_function_workspace_when_debugging()
    switch Input.debug
        case {'yes','y'}
            p = mfilename('fullpath');
            I=strfind(p,'\');
            if nnz(I)==0
                I=strfind(p,'/'); % darn macs
            end
            filename=['debug_' p((I(end)+1):end)];
            save(filename)
        otherwise
    end
end
%Field.Einc =exp(1i*pi/2)*Field.Einc ;
%Field.k = Field.k;
end