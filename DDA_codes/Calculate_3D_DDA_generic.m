 function [P_out,C_out,Field,Array]=Calculate_3D_DDA_generic(Input,Array)
% 3D DDA calculation, using standard formulation for solving
% the self-consistent polarization. Field looped to save memory. mjh-2014
% 1. Make the field E (for given wavelength)
% 2. Make the interaction matrix A
% 3. Solve polarization                        Out: P
% 4. Calculate exctintion etc. cross-sections. Out: C.ext, C.abs etc.
%
% outputs:      P = Polarization matrix
%               C = struct containing the extinction, scattering and
%                   absorption coefficients

Save_function_workspace_when_debugging()

if Input.dimensionality == "cyl_TE"
    error('Calculate_3D_DDA_generic:cylTE', ['cyl_TE is solved by the eigen solver only: set ' ...
        'Input.inversion_tool = ''eig'' (it also returns the driven response P = A\\Einc).']);
end
if contains(string(Input.inversion_tool), "fft") && isfield(Input,'layers') ...
        && isfield(Input.layers,'enable') && Input.layers.enable
    error('Calculate_3D_DDA_generic:fftLayers', ['%s does not support the layered background ' ...
        '(Input.layers): its FFT assumes a homogeneous host. Use eig, gmres, mldivide or linsolve.'], ...
        Input.inversion_tool);
end

%% lighten variables to ease memory issues:
P_tmp = zeros(Input.N,1,'like',Array.r_i);
Array.P_prev=zeros(Input.N,1,'like',Array.r_i);
switch Input.dimensionality
    case '3D'
        P_tmp = zeros(3*Input.N,1,'like',Array.r_i); %  P_tmp = zeros(3*Input.N,1,Array.k_N,'like',Array.r_i);
        P = zeros(3*Input.N,1,Input.range.k_N,'single');
        Array.P_prev=zeros(3*Input.N,1,'like',Array.r_i);
    case '2D'
        P_tmp = zeros(2*Input.N,1,'like',Array.r_i);
        P = zeros(2*Input.N,1,Input.range.k_N,'single');
        Array.P_prev=zeros(2*Input.N,1,'like',Array.r_i);
    case '1Dxcommentedout'
        P_tmp = zeros(Input.N,1,'like',Array.r_i); 
        P = zeros(Input.N,1,Array.k_N,'single');
        Array.P_prev=zeros(Input.N,1,'like',Array.r_i);
end
if nnz(strfind(Input.inversion_tool,'gpu'))>0
    if nnz(strfind(Input.inversion_tool,'fft'))>0
       Array.r_i_conv=gpuArray(Array.r_i_conv);
    end
end
fftw('planner','patient');    % run fftw to speed ffts
if contains(string(Input.inversion_tool), "fft")
    has_mirror = isfield(Input,'mirror_A') && any(Input.mirror_A);
    has_pbc    = isfield(Input,'periodic_boundary') && any(Input.periodic_boundary);
    if has_mirror || has_pbc
        error('Calculate_3D_DDA_generic:fftBoundaries', ['%s does not support mirrors or periodic ' ...
            'boundaries (its FFT assumes a plain finite grid). Use gmres, mldivide, linsolve or eig.'], ...
            Input.inversion_tool);
    end
end

%% make r_ij and R functions for A calculations
Array.r_ij_fun = @(r_i,ii) repmat(r_i(:,ii),[1 size(r_i,1)])-...  % x-component
                           repmat(r_i(:,ii).',[size(r_i,1) 1]);   % y-component
       
Array.R_fun = @(r_i) sqrt(Array.r_ij_fun(r_i,1).^2 + Array.r_ij_fun(r_i,2).^2);   % r_z == 0 identically in XY-plane
t0 = tic;
if Input.inversion_tool == "inv"
    switch Input.dimensionality
        case '3D'
            AinvMaster = single(zeros([prod(Input.arrayN)*3,prod(Input.arrayN)*3,length(Input.range.eV)])); %placeholder
        case '2D'
            AinvMaster = single(zeros([prod(Input.arrayN)*2,prod(Input.arrayN)*2,length(Input.range.eV)]));
        case {'1Dy', '1Dx', '1Dz', 'cyl_TM'}
            AinvMaster = single(zeros([prod(Input.arrayN),prod(Input.arrayN),length(Input.range.eV)]));
    end
end
    
toc
% Calculate NLO DDA for each wavelength kk:
loop_counter=1;modality_timer_tag=' A inverted in ';
Arrayk = Array.k;
if isfield(Array,'r_i_conv'), Arrayr_i_conv = Array.r_i_conv; else, Arrayr_i_conv = []; end
ArrayP_prev = Array.P_prev;
if Input.useparallel
    poolobj = parpool("threads");%"threads"
end
tic 
for kk=1:Input.range.k_N %parfor
    %for parallel computing, P_tmp(:,1) --> P_tmp
   % make N x1 vector too out of the diagonal components:  
   [alpha_diag]=Make_alpha(Array,Input,kk);
   % 1. Define incident field
        [Field]=Make_3D_Field_wl_looped(Input,Array,kk);
        %[Field]=Make_generic_fft_Field_gpu(Input,Array,kk,alpha_diag);

   if kk==1
        disp(['Electric field: [' num2str(Field.Einc(1:3).'), ']. Gaussian: ' num2str(Input.Gaussian.bool)])
   end
   %% Matrix inversion parts:
   switch Input.inversion_tool
       case 'fft_gmres' 
           
           [x0,~]=Make_x0(Field.Einc,ArrayP_prev,kk);                        % initial guess x0
           switch Input.dimensionality
               case '3D'
                   [Afft,siz]=Make_Afft_for_AxGain(kk,Array,Arrayr_i_conv,Input);  % calculate the Afft
                   [P_tmp,flag,Res]=gmres(@(x) GPU_sub_afun(x,Afft,siz,alpha_diag,Input),...
                       Field.Einc(:,1),...
                       Input.gmres.restart,...
                       Input.gmres.tol,...
                       Input.gmres.maxit,...
                       @(x) (1./alpha_diag).\x,[],x0);
               case {'1Dx','1Dy','1Dz'}
                   [Afft,siz]=Make_Afft_for_AxGain(kk,Array,Arrayr_i_conv,Input);  % calculate the Afft
                   [P_tmp,flag,Res]=gmres(@(x) GPU_sub_afun1D(x,Afft,siz,alpha_diag,Input),...
                       Field.Einc(:,1),...
                       Input.gmres.restart,...
                       Input.gmres.tol,...
                       Input.gmres.maxit,...
                       @(x) (1./alpha_diag).\x,[],x0);
               case {"cyl_TM"}
                   [Afft,siz]=Make_cyl_TM_matrix_Afft_SI(kk,Array,Arrayr_i_conv,Input);  % calculate the Afft
                   [P_tmp,flag,Res]=gmres(@(x) GPU_sub_afun1D(x,Afft,siz,alpha_diag,Input),...
                       Field.Einc(:,1),...
                       Input.gmres.restart,...
                       Input.gmres.tol,...
                       Input.gmres.maxit,...
                       @(x) (1./alpha_diag).\x,[],x0);
           end
           P(:,kk) = P_tmp;

       case 'gmres'
           switch  Input.sample.type
               case 'random_surface'     % 'lorentzian Au', 'ring', 'oligomer', 'cm', 'ldr'
                   [Field]=Calculate_Roberts_input_field(Input,Array,P_tmp,kk); 
               otherwise
                   [Field]=Make_3D_Field_wl_looped(Input,Array,kk);
           end
           [x0,~]=Make_x0(Field.Einc,ArrayP_prev,kk);             % initial guess x0
           
           %% 2. Make interaction matrix A for single wavelength
           A = dense_matrix(Input,Array,kk);      % same SI matrix as the eigen solver
           
           %% 3. Solve polarization
           % [P_tmp,flag,Res]=gmres(A,Field.Einc(:,1),...
           %     Input.gmres.restart,...
           %     Input.gmres.tol,...
           %     Input.gmres.maxit,...
           %     @(x) (1./alpha_diag).\x,[],x0);
            [P_tmp,flag,Res]=gmres(A,Field.Einc(:,1),...
               Input.gmres.restart,...
               Input.gmres.tol,...
               Input.gmres.maxit,...
              [],[],x0);
           
           ArrayP_prev=P_tmp;       % store P for next initial quess x0
           P(:,kk)=P_tmp;
           
       case 'fft_gmres_gpu'
           [Afft,siz]=Make_Afft_for_Ax(Arrayk(kk),Arrayr_i_conv,Input);  % calculate the Afft
           [x0,~]=Make_x0(Field.Einc,ArrayP_prev,kk);                        % initial guess x0
                      
           %% try scaling everything closer to unity and then rescaling back
           [P_tmp,flag,Res]=gmres(@(x) gather(GPU_sub_afun(x,Afft,siz,alpha_diag,Input)),...
               gather(Field.Einc(:,1)),...
               Input.gmres.restart,...
               Input.gmres.tol,...
               Input.gmres.maxit,...
               @(x) gather((1./alpha_diag).\x),[],x0);
          
           P(:,kk)=gather(P_tmp);
           ArrayP_prev=P_tmp;       % store P for next initial quess x0

       case 'bicgs_gpu'                      
           [Afft,siz]=Make_Afft_for_Ax(Arrayk(kk),Arrayr_i_conv,Input);  % calculate the Afft
           % Actual iteration line using bicg-stab:
           [P(:,kk),flag,Res]=bicgstab(@(x) gather(GPU_sub_afun(x,Afft,siz,alpha_diag,Input)),...
                    gather(Field.Einc(:,1)),...
                    Input.gmres.tol,...
                    Input.gmres.maxit,...
                    @(x) gather((1./alpha_diag).\x));  
                
       
       case 'gmres_gpu'
           switch  Input.sample.type
               case 'random_surface'     % 'lorentzian Au', 'ring', 'oligomer', 'cm', 'ldr'
                   [alpha_diag]=Make_alpha(Array,Input,1);
                   if length(Input.arrayN) >2    
                       [Field]=Calculate_3D_robert_input_field(Input,Array,kk);
                   else
                       [Field]=Calculate_Roberts_input_field(Input,Array,P_tmp,kk);
                   end
                  
               otherwise
                   [Field]=Make_3D_Field_wl_looped(Input,Array,kk);
           end           
           [x0,~]=Make_x0(Field.Einc,ArrayP_prev,kk);             % initial guess x0
           
           %% 2. Make interaction matrix A for single wavelength
           A = dense_matrix(Input,Array,kk);      % same SI matrix as the eigen solver

           %% 3. Solve polarization
           [P_tmp,flag,Res]=gmres(A,Field.Einc(:,1),...
               Input.gmres.restart,...
               Input.gmres.tol,...
               Input.gmres.maxit,...
               @(x) gather((1./alpha_diag).\x),[],x0);
                      
           ArrayP_prev=P_tmp;       % store P for next initial quess x0
           P(:,kk)=gather(P_tmp);
       case 'inv'
            %Inverting tool = inv
            % Just inverse the A matrix. Probably extremly slow
    
            [Field]=Make_3D_Field_wl_looped(Input,Array,kk);
                A = dense_matrix(Input,Array,kk);      % same SI matrix as the eigen solver
                %{
                disp("Check if hermitian")
                if ishermitian(A)
                     disp("The matrix A is hermitian.")
                else
                    disp("The matrix is not hermitian")
                end
                %}
                
                %Ainv = inv(A)
                Ainv = A\eye(size(A));
                P_tmp = Ainv*Field.Einc(:,1);
                P(:,kk)=P_tmp;
                [~]=Display_calc_times_and_convergence_specs('fwm',kk,Array.k,modality_timer_tag);
                flag = 0; Res = 0;
                AinvMaster(:,:,kk) = Ainv;
       case {'mldivide','linsolve'}
            % Direct dense solve of A*P = Einc
            [Field]=Make_3D_Field_wl_looped(Input,Array,kk);
            A = dense_matrix(Input,Array,kk);
            if Input.inversion_tool == "linsolve"
                P_tmp = linsolve(A, double(Field.Einc(:,1)));
            else
                P_tmp = A \ double(Field.Einc(:,1));
            end
            P(:,kk) = P_tmp;
            flag = 0; Res = 0;
       otherwise
           error('Calculate_3D_DDA_generic:tool', 'Unknown Input.inversion_tool "%s".', Input.inversion_tool)
   end
   
   % Counter to check if result converged, if there was a problem, flag's written into counter and error displayed:
   [~]=Display_calc_times_and_convergence_specs('convergence',kk,Arrayk,[' A(w) inverted in '],flag,Res);        
   %clear x
   Display_calc_times_and_convergence_specs('fwm',kk,Array.k,modality_timer_tag, 0,NaN(1), t0);
end 
if Input.useparallel
    delete(poolobj);
end
toc
eps0 =  8.8541878188e-12; %F/m
[Field]=Make_3D_Field_wl_looped(Input,Array,kk);

for kk = 1:Input.range.k_N
    [alpha_diag]=Make_alpha(Array,Input,kk);
    [Field]=Make_3D_Field_wl_looped(Input,Array,kk);
   
   %% 4. Extinction and absorption, both per array area and incident intensity
   % (dda_extinction has the same normalization). For SI polarizabilities the
   % prefactor is 1/(4*pi*eps0*n^2) and the radiation reaction k^3/(6*pi*eps0*n^2)
   % (k^2/(4*eps0*n^2) for cylinders); for Gaussian ones 1 and 2/3*k^3.
   switch string(Input.dimensionality)
       case "3D",  nc = 3;
       case "2D",  nc = 2;
       otherwise,  nc = 1;
   end
   Ep = reshape(double(Field.Einc(:,1)), nc, []);
   I0 = mean(sum(abs(Ep).^2, 1));                         % per particle
   if Input.dimensionality == "cyl_TM", cell_area = Input.p(1);
   else,                                cell_area = Input.p(1)*Input.p(2);
   end
   kk_k = Array.k(kk);  nn = Array.n(min(kk, numel(Array.n)));
   pref = 4*pi*kk_k/(Input.N*cell_area*I0);
   if dda_alpha_is_SI(Input,Array)
       pref = pref/(4*pi*eps0*nn^2);
       if Input.dimensionality == "cyl_TM", rad = kk_k^2/(4*eps0*nn^2);
       else,                                rad = kk_k^3/(6*pi*eps0*nn^2);
       end
   else
       rad = 2/3*kk_k^3;
   end
   Pk_ = double(P(:,kk));  Ek_ = double(Field.Einc(:,1));  ad_ = double(alpha_diag(:));
   ext(kk)    = pref*sum(imag(Pk_.*conj(Ek_)));
   absorp(kk) = pref*sum(imag(Pk_.*conj(1./ad_).*conj(Pk_)) - rad*abs(Pk_).^2);
   for xx = 1:nc
       ext_vec(xx,kk)     = pref*sum(imag(Pk_(xx:nc:end).*conj(Ek_(xx:nc:end))));
       abs_vec(xx,kk)     = pref*sum(imag(Pk_(xx:nc:end).*conj(1./ad_(xx:nc:end)).*conj(Pk_(xx:nc:end))) ...
                                     - rad*abs(Pk_(xx:nc:end)).^2);
       ext_vec_max(xx,kk) = pref*Input.N*max(imag(Pk_(xx:nc:end).*conj(Ek_(xx:nc:end))));
       abs_vec_max(xx,kk) = pref*Input.N*max(imag(Pk_(xx:nc:end).*conj(1./ad_(xx:nc:end)).*conj(Pk_(xx:nc:end))) ...
                                     - rad*abs(Pk_(xx:nc:end)).^2);
   end
   %% store also the local field, or the mean of it etc.
   %
   if Input.Eloc
       
       switch Input.dimensionality
           %
           case {'1Dx','1Dy','1Dz','cyl_TM'}
               E_enhancement_sum(1,kk) = squeeze(sum(sum(reshape( P(:,kk)./(Field.Einc(:,1).*alpha_diag(:)),Input.arrayN(1),Input.arrayN(2)))))/(Input.arrayN(1)*Input.arrayN(2));
               E_enhancement_max(1,kk) = squeeze(max(max(reshape( P(:,kk)./(Field.Einc(:,1).*alpha_diag(:)),Input.arrayN(1),Input.arrayN(2))))); % /(Input.arrayN^2); 
           %}
           case '3D'
               for ii=1:3
                   if length(Input.arrayN) < 2
                       E_enhancement_sum(ii,kk) = squeeze(sum(sum(sum(reshape( P(ii:3:end,kk)./(Field.Einc(ii:3:end,1).*alpha_diag(ii:3:end)),Input.arrayN,Input.arrayN)))))/(Input.N);
                       E_enhancement_max(ii,kk) = squeeze(max(max(reshape( P(ii:3:end,kk)./(Field.Einc(ii:3:end,1).*alpha_diag(ii:3:end)),Input.arrayN,Input.arrayN))));                       
                   else
                       E_enhancement_sum(ii,kk) = squeeze(sum(sum(sum(reshape( P(ii:3:end,kk)./(Field.Einc(ii:3:end,1).*alpha_diag(ii:3:end)),Input.arrayN(1),Input.arrayN(2),Input.arrayN(3))))))/(Input.N);
                       E_enhancement_max(ii,kk) = squeeze(max(P(ii:3:end,kk)./(Field.Einc(ii:3:end,1).*alpha_diag(ii:3:end)))); % /(Input.N)
                   end
               end
       end
   end
   %}
   % 5. display calculation times
   
   
end
C.ext = ext; C.abs=absorp;C.ext_vec=ext_vec; C.abs_vec = abs_vec;
C.solver_version = 4;   % as Calculate_DDA_eigenmodes: 4 = mirror fix in Make_int_matrix_A_SI
C.ext_vec_max = ext_vec_max; C.abs_vec_max = abs_vec_max;

         
            

C.sca = C.ext-C.abs;    % scattering coefficient
C.sca_vec = C.ext_vec-C.abs_vec;    % scattering coefficient
C.sca_vec_max = C.ext_vec_max-C.abs_vec_max;    % scattering coefficient

%% gather and clear variables for storing
% gather variables out from GPU if needed:
if nnz(strfind(Input.inversion_tool,'gpu'))>0
    C_out.sca = gather(C.sca);
    C_out.abs = gather(C.abs);
    C_out.ext = gather(C.ext); 
    C_out.sca_vec = gather(C.sca_vec);
    C_out.abs_vec = gather(C.abs_vec);
    C_out.ext_vec = gather(C.ext_vec);
    
    C_out.sca_vec_max = gather(C.sca_vec_max);
    C_out.abs_vec_max = gather(C.abs_vec_max);
    C_out.ext_vec_max = gather(C.ext_vec_max);
       
    P_out=gather(P);
    
    if Input.Eloc
        C_out.E_loc_sum=gather(E_enhancement_sum);
        C_out.E_loc_max=gather(E_enhancement_max);
    end
else
    C_out.sca = C.sca;
    C_out.abs = C.abs;
    C_out.ext = C.ext;
    C_out.sca_vec = C.sca_vec;
    C_out.abs_vec = C.abs_vec;
    C_out.ext_vec = C.ext_vec;
    
    C_out.sca_vec_max = C.sca_vec_max;
    C_out.abs_vec_max = C.abs_vec_max;
    C_out.ext_vec_max = C.ext_vec_max;
    
    P_out=P;
   
    if Input.Eloc
        C_out.E_loc_sum=E_enhancement_sum;
        C_out.E_loc_max=gather(E_enhancement_max);
    end
end
  
if find(isfield(Array,'alpha'))==1
    Array.alpha1=squeeze(Array.alpha(1,:,:,:));
end

% remove 'r_ij' and 'alpha' 'r_i_conv' fields:
fields={'r_ij','alpha','r_i_conv','P_prev','r_i'}; 
fields={'r_ij','alpha','r_i_conv','P_prev'}; 

I=find(isfield(Array,fields));
if nnz(I)
    for ii=I
        Array=rmfield(Array,fields(ii));
    end
end

%%
% Function to save the function workspace when debugging
    function Save_function_workspace_when_debugging
        if Input.debug
            p = mfilename('fullpath');
            I=strfind(p,'\');
            if nnz(I)==0
                I=strfind(p,'/'); % darn macs
            end
            filename=['debug_' p((I(end)+1):end)];
            save(filename)
        end
    end


end

function [x0,scaler]=Make_x0(Einc,P_prev,kk)
% Makes an educated initial guess, x0(kk=1) = Zeros, x0(kk~=1) = P_prev.
% Also finds a scaler to normalize parameters, to reduce numerical error.
if kk==1
    x0=zeros(size(P_prev),'like',P_prev); % zero vector to begin with
else
    x0=P_prev;        % Make an educated initial guess
end
scaler=max(Einc);     % scale variables to reduce numerical error in iteration:
end

function [Afft,siz]=Make_Afft_for_Ax(k,r_i,Input)
% Calculate convolution form of Aij for single wavelength

A_BT=cell(3,3);         % BT = Block Toeplitz 
I3=eye(3,'like',r_i);   % construct eye matrix into same memory as r_i variable

%% Define variables to speed-up calculations of A_BT:
R=sqrt(r_i(:,1).^2 + r_i(:,2).^2 + r_i(:,3).^2);
R(R<1e-14)=0;   % force zero dipole R to zero, so that isinf/isnan do find the right spot after the for loop

exp_term=exp(1i*k.*R)./(R);         
R_sqrd_term=(1i*k.*R-1)./(R.^2);

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

          A_Arr = exp_term.*(k^2.*(r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)) ...
                  +R_sqrd_term.*(3.*r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)));  
              
          % find NaN location, as that's where the main dipole is and input 1/alpha into it:
          A_Arr(isnan(A_Arr))=1e-10;  % 1/1e-30 was 1e-30 before 10.3.2017
          A_Arr(isinf(A_Arr))=1e-10;  % 1/1e-30
          
          A_BT{ii,jj}=reshape(A_Arr,[Input.Nmax(1)*2-1 Input.Nmax(2)*2-1  Input.Nmax(3)*2-1]);
    end            
end
clear A_Arr exp_term R_sqrd_term r_i_per_R

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
        Afft{ii,jj}=fftn(A_BT{ii,jj},siz.total);
    otherwise
        for jj=ind_jj
            for ii=ind_ii
                Afft{ii,jj}=fftn(A_BT{ii,jj},siz.total);     % (jj,z)
            end
        end
end

end

function [y]=GPU_sub_afun(x,Afft,siz,alpha,Input)
%% Performs fast matrix-vector multiplication related to iteratively solving
% source polarization in DDA:  
%             Ax, where A = interaction matrix 
%                       x = trial of source polarization
%
% Script is based on the fact that A is Block Toeplitz (BT) for rectangular
% dipole arrays, and hence the multiplication can be turned into a
% convolution. The convolution is calculated using the regular fft-based
% speed-up schemes. More info can be found from e.g. 
%
% Draine papers, Flatau Review on DDA 
%
% In addition, only small part of the A has to computed and
% convolved as there's a lot of redundancy in the full A matrix for the
% system. 
%
% correct size for the interaction matrix is 2N-1, if the traditional A is NxN
% But we only need to calculate a part of A, where the dipole has N-1 neighbours in 
% all directions. 
% e.g. N=10, and NxN array traditionally requires constructing A of size N^2xN^2 
%
% Things to make into the codes, 
% 1) check that indeed we're working with square/rectangular arrays
% Input.arrayN  or Input.N
%
% How the actual algorithm works:
% 1) construct the right A_block for convolution
% 2) check sizes and make ffts
% 3) make ifft and crop the right part of y
%
% 1) For NxN array, we need to calculate A_block = B, which corresponds to 
%    an interaction matrix A of a dipole surrounded by N-1 neighbouring dipoles
%    in all directions.
%
% x=gpuArray(single(x));
if nnz(strfind(Input.inversion_tool,'gpu'))>0
    x=gpuArray(single(x));
end

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
X_BT=cell(length(ind_ii),1);

%% 2.  structure x:
for ii=ind_ii  
    % Zero padding and reshaping of x in to X 
    Ntot=Input.Nmax(1)*Input.Nmax(2)*Input.Nmax(3);    
    switch Input.dimensionality
        case {'1Dx','1Dy','1Dz'}
            X_BT{ii}=single(reshape(x(1:Ntot),Input.Nmax(1),Input.Nmax(2),Input.Nmax(3)));
        otherwise
            X_BT{ii}=single(reshape(x(ii:3:3*Ntot),Input.Nmax(1),Input.Nmax(2),Input.Nmax(3)));
    end           
end
 
%% 3. Fourier and inverse fourier transforms to perform the convolution:
% OBS! fftn might be able to do an even higher dimensionality FFT faster than with now is being for looped! 10.12.2015 mjh
y=zeros(size(x),'like',alpha);
y_full=cell(length(ind_jj),1);
Xfft = cell(length(ind_jj),1);
for jj=ind_jj
    Xfft{jj}=fftn(X_BT{jj},siz.total);  % (jj,x)
end

%% 4. Calculate matrix Aij for single wavelength
switch Input.dimensionality
    case {'1Dx','1Dy','1Dz'}                        
        %y_full=ifftn(Afft{ii,jj}.*fftn(X_BT{ii},siz.total));          % actual convolution part
        y_full=ifftn(Afft{ii,jj}.*Xfft{ii});          % actual convolution part
        y=reshape(y_full(ceil((siz.A(1)-1)/2)+(1:siz.X(1)),ceil((siz.A(2)-1)/2)+(1:siz.X(2)),:),[],1);
    otherwise
        for jj=ind_jj
            % actual convolution part
            y_full{jj}=ifftn(Afft{1,jj}.*Xfft{1}+...   % (jj,x)
                             Afft{2,jj}.*Xfft{2}+...   % (jj,y)
                             Afft{3,jj}.*Xfft{3});     % (jj,z)
            
            % part of y, the A_ij * P_j part, taken as properly reshaped center-part of the convolution:
            y(jj:3:end)=reshape(y_full{jj}(ceil((siz.A(1)-1)/2)+(1:siz.X(1)),ceil((siz.A(2)-1)/2)+(1:siz.X(2)),ceil((siz.A(3)-1)/2)+(1:siz.X(3))),[],1);            
        end
end

%% 5. Form y = alpha_i^-1 * P_i   +   A_ij * P_j , which is the full trial
y = (1./alpha).*x + y;

end



function [A]=Make_1D_int_matrix_A_sub_varying_gain(Input,Array,kk)
% Form A = N x N block interaction matrix, where i=j terms are inverses of
% single particle polarizabilities alpha
    
%% For calculating R
    function R=R_fun(r_i)
        R=sqrt(Array.r_ij_fun(r_i,1).^2 + Array.r_ij_fun(r_i,2).^2);   % r_z == 0 identically in XY-plane
        
        % replace zero R terms with large distance so, that those won't affect results             
        R(1:size(R,1)+1:end) = 1e6;     
    end
%% For calculating R_hat
    function R_jk_hat=R_jk_hat_fun(r_i,ii)
        if ii<3
            R_jk_hat =  Array.r_ij_fun(r_i,ii)./R_fun(r_i);
        else
            R_jk_hat=zeros(size(r_i,1));      % r_z == 0 identically in XY-plane
        end
    end
%% 
r_i = Array.r_i;
R=sqrt(r_i(:,1).^2 + r_i(:,2).^2 + r_i(:,3).^2);
R(R<1e-14)=0;   % force zero dipole R to zero, so that isinf/isnan do find the right spot after the for loop
r_i_per_R(:,1)=r_i(:,1)./R;
r_i_per_R(:,2)=r_i(:,2)./R;
r_i_per_R(:,3)=r_i(:,3)./R;

%% Calculate off-diagonal blocks for matrix Aij for single wavelength
k=Array.k(kk);

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
end
if exist('Input.gain')
    gain_term = exp(Input.gain*k*R);
else
    gain_term = 1;
end

exp_term=exp(1i*k.*R)./(R);         
R_sqrd_term=(1i*k.*R-1)./(R.^2);
% Form the A and Aself matrices:
% A=zeros(size(R_fun(Array.r_i),1),size(R_fun(Array.r_i),1));
for ii=ind_ii
    for jj=ind_jj
            % 26.1.2017 was like this before.. problem with the ii,jj
            % indices? -mjh s
%           A_Arr = exp_term.*(k^2.*(r_i_per_R(:,ii).*r_i_per_R(:,ii)-I3(ii,jj)) ...
%                   +R_sqrd_term.*(3.*r_i_per_R(:,ii).*r_i_per_R(:,ii)-I3(ii,jj)));  

          A = gain_term.*exp_term.*(k^2.*(r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)) ...
                  +R_sqrd_term.*(3.*r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)));  
              
          % find NaN location, as that's where the main dipole is and input 1/alpha into it:
          A(isnan(A))=1e-10;  % 1/1e-30 was 1e-30 before 10.3.2017
          A(isinf(A))=1e-10;  % 1/1e-30
    end            
end
A =reshape(A,[Input.Nmax(1) Input.Nmax(2)  Input.Nmax(3)]);
%{
if (ind_ii==3 && ind_jj==3)
    %% A only function of R etc.
    A = exp(1i*k.*R_fun(Array.r_i))./(4*pi*R_fun(Array.r_i)).*(k^2.*(-1)+...
        (1i*k.*R_fun(Array.r_i)-1)./(R_fun(Array.r_i).^2).*(-1))+...
        diag(repmat(squeeze(1./Array.alpha(1,3,3,kk)),Input.N,1));       % A_ii , i.e. polarizabilities
else
    %% R_hats are needed in calculation of A
    for ii=ind_ii     % loop over the right 3 x 3 indices of the dyadic
        for jj=ind_jj
            A = exp(1i*k.*R_fun(Array.r_i))./(4*pi*R_fun(Array.r_i)).* ...
                (k^2.*(R_jk_hat_fun(Array.r_i,ii).*R_jk_hat_fun(Array.r_i,jj)-I3(ii,jj)) ...
                +(1i*k.*R_fun(Array.r_i)-1)./(R_fun(Array.r_i).^2).*...
                (3*R_jk_hat_fun(Array.r_i,ii).*R_jk_hat_fun(Array.r_i,jj))-I3(ii,jj))+...
                diag(repmat(squeeze(1./Array.alpha(1,ind_ii,ind_jj,kk)),Input.N,1));
        end
    end
end

nnzNaN=find(isnan(A));
A((nnzNaN))=0;
%}
end

function [Afft,siz]=Make_Afft_for_Ax_0603_2016_backup(k,r_i,Input)
% Calculate matrix Aij for single wavelength

I3=eye(3,'like',r_i);   % construct eye matrix into same memory as r_i variable
A_BT=cell(3,3);

R=sqrt(r_i(:,1).^2 + r_i(:,2).^2 + r_i(:,3).^2);
R(R<1e-14)=0;   % force the zero dipole R to be zero, so that isinf/isnan do find the right spot after the for loop
exp_term=exp(1i*k.*R)./(R);        % trial to speed up calculation of A 
R_sqrd_term=(1i*k.*R-1)./(R.^2);
r_i_per_R(:,1)=r_i(:,1)./R;
r_i_per_R(:,2)=r_i(:,2)./R;
r_i_per_R(:,3)=r_i(:,3)./R;

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

%% 5.4.2016 consideration, should the r_hat * r_hat be rather a 3x3 array and not just a diag(3x3) ?
% % % 1. loop over the right 3 x 3 indices of the dyadic to calculate array A
% % for ii=ind_ii
% %     for jj=ind_jj
% %           A_Arr = exp_term.*(k^2.*(r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)) ...
% %                   +R_sqrd_term.*(3.*r_i_per_R(:,ii).*r_i_per_R(:,jj)-I3(ii,jj)));  
% %           
% %           % find NaN location, as that's where the main dipole is and input 1/alpha into it:
% %           A_Arr(isnan(A_Arr) | isinf(A_Arr))=1e-30;  % 1/1e-30          
% %           A_BT{ii,jj}=reshape(A_Arr,[Input.Nmax(1)*2-1 Input.Nmax(2)*2-1  Input.Nmax(3)*2-1]);
% %     end            
% % end

% 1. loop over the right 3 x 3 indices of the dyadic to calculate array A
for ii=ind_ii
    for jj=ind_jj
          A_Arr = exp_term.*(k^2.*(r_i_per_R(:,ii).*r_i_per_R(:,ii)-I3(ii,jj)) ...
                  +R_sqrd_term.*(3.*r_i_per_R(:,ii).*r_i_per_R(:,ii)-I3(ii,jj)));  
          
          % find NaN location, as that's where the main dipole is and input 1/alpha into it:
          A_Arr(isnan(A_Arr) | isinf(A_Arr))=1e-30;  % 1/1e-30          
          A_BT{ii,jj}=reshape(A_Arr,[Input.Nmax(1)*2-1 Input.Nmax(2)*2-1  Input.Nmax(3)*2-1]);
    end            
end
clear A_Arr exp_term R_sqrd_term

% Check the sizes for making the correct ffts
siz.A = [size(A_BT{ii}) 1];
siz.X = [Input.Nmax(1) Input.Nmax(2) Input.Nmax(3)];
siz.total =  2.^nextpow2(siz.A + siz.X -1);
% siz.total =  siz.A + siz.X -1;

%% Fourier and inverse fourier transforms to perform the convolution:
%% OBS! fftn might be able to do an even higher dimensionality FFT faster than with now is being for looped! 10.12.2015 mjh
Afft=cell(length(ind_ii),length(ind_jj));
switch Input.dimensionality
    case {'1Dx','1Dy','1Dz'}
        Afft{ii,jj}=fftn(A_BT{ii,jj},siz.total);
    otherwise
        for jj=ind_jj
            for ii=ind_ii
                Afft{ii,jj}=fftn(A_BT{ii,jj},siz.total);     % (jj,z)
            end
        end
end

end

function [y]=afun_Ax_fft_0603_2016_backup(x,Afft,siz,alpha,Input)
%% Performs fast matrix-vector multiplication related to iteratively solving
% source polarization in DDA:  
%             Ax, where A = interaction matrix 
%                       x = trial of source polarization
%
% Script is based on the fact that A is Block Toeplitz (BT) for rectangular
% dipole arrays, and hence the multiplication can be turned into a
% convolution. The convolution is calculated using the regular fft-based
% speed-up schemes. More info can be found from e.g. 
%
% Draine papers, Flatau Review on DDA 
%
% In addition, only small part of the A has to computed and
% convolved as there's a lot of redundancy in the full A matrix for the
% system. 
%
%
% correct size for the interaction matrix is 2N-1, if the traditional A is NxN
% But we only need to calculate a part of A, where the dipole has N-1 neighbours in 
% all directions. 
% e.g. N=10, and NxN array traditionally requires constructing A of size N^2xN^2 
%
% Things to make into the codes, 
% 1) check that indeed we're working with square/rectangular arrays
% Input.arrayN  or Input.N
%
% How the actual algorithm works:
% 1) construct the right A_block for convolution
% 2) check sizes and make ffts
% 3) make ifft and crop the right part of y
%
% 1) For NxN array, we need to calculate A_block = B, which corresponds to 
%    an interaction matrix A of a dipole surrounded by N-1 neighbouring dipoles
%    in all directions.
%

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
X_BT=cell(length(ind_ii),1);

% 2.  structure x:
for ii=ind_ii  
    %% Zero padding and reshaping of x in to X 
    Ntot=Input.Nmax(1)*Input.Nmax(2)*Input.Nmax(3);
    
    switch Input.dimensionality
        case {'1Dx','1Dy','1Dz'}
            X_BT{ii}=single(reshape(x(1:Ntot),Input.Nmax(1),Input.Nmax(2),Input.Nmax(3)));
        otherwise
            X_BT{ii}=single(reshape(x(ii:3:3*Ntot),Input.Nmax(1),Input.Nmax(2),Input.Nmax(3)));
    end
            
end
 
%% Fourier and inverse fourier transforms to perform the convolution:
%% OBS! fftn might be able to do an even higher dimensionality FFT faster than with now is being for looped! 10.12.2015 mjh
y=zeros(size(x),'like',alpha);
y_full=cell(length(ind_jj),1);
Xfft = cell(length(ind_jj),1);

for jj=ind_jj
    Xfft{jj}=fftn(X_BT{jj},siz.total);  % (jj,x)
end

%% Calculate matrix Aij for single wavelength
switch Input.dimensionality
    case {'1Dx','1Dy','1Dz'}                        
        y_full=ifftn(Afft{ii,jj}.*fftn(X_BT{ii},siz.total));          % actual convolution part
        y=reshape(y_full(ceil((siz.A(1)-1)/2)+(1:siz.X(1)),ceil((siz.A(2)-1)/2)+(1:siz.X(2)),:),[],1);
    otherwise
        for jj=ind_jj
            % actual convolution part
            y_full{jj}=ifftn(Afft{1,jj}.*Xfft{1}+...   % (jj,x)
                             Afft{2,jj}.*Xfft{2}+...   % (jj,y)
                             Afft{3,jj}.*Xfft{3});     % (jj,z)
            
            %% part of y, the A_ij * P_j part, taken as properly reshaped center-part of the convolution:
            y(jj:3:end)=reshape(y_full{jj}(ceil((siz.A(1)-1)/2)+(1:siz.X(1)),ceil((siz.A(2)-1)/2)+(1:siz.X(2)),ceil((siz.A(3)-1)/2)+(1:siz.X(3))),[],1);            
        end
end

%% y = alpha_i^-1 * P_i   +   A_ij * P_j , which is the full trial
y = (1./alpha).*x + y;

end

%% ------------------------------------------------------------------------
function A = dense_matrix(Input, Array, kk)
% The dense interaction matrix, built exactly as for the eigen solver: SI
% units, 1D/2D/3D/cyl_TM, mirrors and periodic boundaries.
if Input.dimensionality == "cyl_TM"
    A = Make_cyl_TM_matrix_A_SI(Input, Array, kk);
else
    A = Make_int_matrix_A_SI(Input, Array, kk);
end
end
