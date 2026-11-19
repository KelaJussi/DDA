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

%% Field profile (Input.preset_field)
%   "none"      plane wave (above)
%   "Gaussian"  Hermite-Gaussian beam HG(mx,my) on top of it (gaussian_beam below)
switch lower(char(string(Input.preset_field)))
    case 'none'
    case {'gaussian', 'gaussian_proper'}            % gaussian_proper: the old name
        Field.Einc = gaussian_beam(Field.Einc, Input, Array_in, jj);
    otherwise
        error('Make_3D_Field_wl_looped:preset', ...
            'Input.preset_field = "%s" is not available; use "none" or "Gaussian".', ...
            char(string(Input.preset_field)));
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

%% ------------------------------------------------------------------------
function E = gaussian_beam(E, Input, Array, jj)
% Hermite-Gaussian beam HG(mx,my), mx = Input.Hermite_order_x, my = Input.Hermite_order_y,
% multiplying the plane wave E at the particle positions:
%   u(x,y) = (w0/w) H_mx(sqrt(2) x/w) H_my(sqrt(2) y/w) exp(-(x^2+y^2)/w^2 - i k (x^2+y^2)/(2R))
%   w = w0 sqrt(1+(z/zR)^2),  1/R = z/(z^2+zR^2),  zR = k w0^2/2 = pi w0^2 n/lambda0
% w0 = Input.Gaussian_w0 (waist), z = Input.Gaussian_z: distance from the array to the
% waist (z > 0: the waist lies beyond the array, the beam converges at the array;
% z = 0: waist on the array, flat phase). The constant Gouy phase is left out.
% x, y are measured from the centre of the array (midway between the outermost
% particles), at the actual particle positions (shifted lattices, edited layouts).
% Cylinders (cyl_TM, cyl_TE, infinite along y): a 1D beam in x, amplitude
% sqrt(w0/w); Hermite_order_y does not apply.
w0 = Input.Gaussian_w0;  z = Input.Gaussian_z;
mx = round(Input.Hermite_order_x);  my = round(Input.Hermite_order_y);
k  = Array.k(jj);                                  % wavenumber in the host medium
zR = k*w0^2/2;
w  = w0*sqrt(1 + (z/zR)^2);
invR = z/(z^2 + zR^2);                             % 0 at the waist
r  = double(Array.r_i);
x  = r(:,1) - (min(r(:,1)) + max(r(:,1)))/2;
y  = r(:,2) - (min(r(:,2)) + max(r(:,2)))/2;
if any(strcmp(char(Input.dimensionality), {'cyl_TM', 'cyl_TE'}))
    u = sqrt(w0/w) * hermite(mx, sqrt(2)*x/w) .* exp(-x.^2/w^2 - 1i*k*x.^2*invR/2);
else
    rho2 = x.^2 + y.^2;
    u = (w0/w) * hermite(mx, sqrt(2)*x/w) .* hermite(my, sqrt(2)*y/w) ...
        .* exp(-rho2/w^2 - 1i*k*rho2*invR/2);
end
ncomp = numel(E)/numel(u);                         % components interleaved per particle
E = reshape(reshape(E, ncomp, []) .* reshape(u, 1, []), [], 1);
end

function h = hermite(n, t)
% physicists' Hermite polynomial H_n(t): H0 = 1, H1 = 2t, H_{m+1} = 2t H_m - 2m H_{m-1}
h = ones(size(t));
if n == 0, return; end
hm = h;  h = 2*t;
for m = 1:n-1
    hn = 2*t.*h - 2*m*hm;  hm = h;  h = hn;
end
end
