function []=BroadBand_NLO_SLRs_pump1000nm_tunigroup()
% iteratively fits the FDTD extinction with the calculated DDA. 
% searches the best scattering cross-secion by comparing results calculated
% using (1) FDTD and (2) DDA.

%% by breaking a short pulse into plane waves incident at different angles,
% SLR-system could be made to enhance responses via narrow resonances that still 
% effectively enhance the signals broadband. People might think this neat.


%% predefined input arguments
Make_predefined_Inputs(); 

Input.inversion_tool='fft_gmres';      % 'fft_gmres_gpu' linsolve, mldivide , gmres_gpu
Input.sample.n=1.51;  
 
%Input.arrayN=[50 50 1];                         % array size is NxN    2
Input.arrayN=[100 100 1];                         % array size is NxN    2
Input.arrayN=[100 100 1];                         % array size is NxN    2
%% anisotropic scattering coefficient:
Input.particle.particle_type='Au_L'; % 'Au_prism'
% Input.particle.scatter_strength=[0.3 0.05 0.1]*1e-5/(4*pi);

%% try increasing scattering cross section
%% with over 5 scatter strength values, MLWA approximation is used.
% [x-pol x-pol y-pol y-pol z-pol z-pol]???
%
Input.particle.scatter_strength=[0.002 0.6 0.9 0.9 0.1 0.1]*1e-5/(4*pi); % gives T=0.79 for y-pol
Input.particle.scatter_strength=[0.002 0.6 1 1 0.1 0.1]*1e-5/(4*pi);

%% MLWA complicates finding fits, since now both the single particle responses and arrays depend on scattering strength..

Input.particle.tau=[0.8 1.15 1.15]*1e-15;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15
Input.particle.tau=[2 0.2 2]*1e-15;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15
% Input.particle.tau=[4 4 4]*1e-15;   % used to be 0.3e-15 for silver, check how could change things 0.05e-15
% Input.particle.lambda_peak_vec=1e-9*[657-200 477-125 1000]; % plasmon resonances lambda_x, lambda_y, lambda_z
Input.particle.lambda_peak_vec=1e-9*[500 330 1000]; % plasmon resonances lambda_x, lambda_y, lambda_z




Py=(710-14)*1e-9;
Px=(662)*1e-9; % Center SLR to 1000 nm
%Py=(695)*1e-9; % Center SLR to 1050 nm

% Px=Py/2+5e-9;        % to make a doubly-resonant array
Py=(400)*1e-9;       % 394.3 resulted in SLR at 596.7 nm to make a doubly-resonant array
Pz=Px; %% WHY?         % 
Input.Pz=Pz;


%Input.arrayN=[round(300e-6/Px) round(300e-6/Py) 1];  %Array corresind to 300 um x300 um array
%% fdtd 0.3943  0792 microns

Input.p(1,1,:)=[Px Py Input.Pz]';

% Changing parameter in calculations is now the period of the array (in x-direction)  
Input.change.paratype='p';           
for ii=1:size(Px,2)
    dip(ii,:,:)=[Px(ii) Py(ii) Pz];
end
Input.change.values=dip;  Input.Px=Px;  Input.Py=Py;

% Plotting handles:
plotting.oneD.visible=true;     % show single 1D spectra during the calculations
plotting.twoD.visible=false;
plotting.limits_2D=[1040 1260]*1e-9; % [980 1320]*1e-9; % these limits are for the fundamental fields.
plotting.limits_2D=[1006 1008]*1e-9; % [980 1320]*1e-9; % these limits are for the fundamental fields.
plotting.type = '1-ext'; %'ext' '1-ext'
Input.plotting=plotting;

fine_wl_range=[945 1105];
fine_wl_reso=1;
%fine_wl_reso=10;
Input.range.wl_tmp=rot90([500:1:fine_wl_range(1) ...
                (fine_wl_range(1)+1):fine_wl_reso:fine_wl_range(2) ...
                (fine_wl_range(2)+1):5:1200]*1e-9,2);
Input.range.wl_tmp = rot90([960:0.4:1040]*1e-9,2);
Input.range.eV=nm2ev(Input.range.wl_tmp);

%% the Ey & kx case %   
Input.theta=1;       % facilitate dispersion calculations, 32 for 2 waves !! refractive index!!
Input.phi=0;         % angle in xy-plane, 0 = along x, rot towards y-axis
Input.psi=0;        % rotation of x'y'-plane, 0: x=x', rot towards y'-axis 

%% Dispersion relation calculation:
%Input.theta_var=0:0.5:5;       % facilitate dispersion calculations, 32 for 2 waves !! refractive index!!
% Input.theta_var=[0:0.25:2];
% Input.theta_var=[0 2 3 4];
%Input.theta_var=[0:0.25:4];
Input.theta_var=0:0.1:1;
%Input.theta_var=0:1:4;
%Input.theta_var=[0 4 5 6 7 8];

for ii=Input.theta_var
    Input.theta=ii;
    
    % make a smarter wl input vector:
    %     DOs occur at 
    %     Input.range.wl_tmp=rot90([800:4.5:1099 1100:0.5:1199 1200:4.5:1300]*1e-9,2);
    NLO_DDA_Dispersion_period_varies(Input,plotting)
end
 

%% NESTED FUNCTIONS %%%%

% PRE-DEFINED inputs, use SI-units  
    function Make_predefined_Inputs        
        warning('off')
        
        Input.project='DDA';
        
        Input.sample.type='lorentzian Au';     % 'lorentzian Au', 'ring', 'oligomer', 'cm', 'ldr'               
        Input.inversion_tool='fft_gmres';      % linsolve, mldivide , gmres_gpu
     
        Input.Eloc=true;
        Input.debug=false;
        Input.dimensionality='3D';            % 1Dz, 3D, 1Dx
        Input.array_type='square';      % lattice type :    'square' 'hexagonal'

        plotting.type = 'ext';           % ext = extinction ,  transmission = 1-ext
        plotting.yaxis = 'lambda';  	 % eV , lambda
        plotting.xaxis = 'p_nm';         % '1/p_nm' 'p_nm'

        % Make rootfolder for the calculations by taking the current folder as the master:
        p = mfilename('fullpath');
        [pathstr,name,~]=fileparts(p);
        Input.code_rootfolder=pathstr;
        Input.m_name=name;
        Input.stl_folder=[Input.code_rootfolder 'stl_samples\'];
            
        % make host computer specific rootfolder for saving data:
        %[Input.sys] = Get_cpu_list_and_folder_hierarchy();
        [Input.sys] = Get_cpu_list_and_folder_hierarchy_TS();
        Input.date=date;
        
        % switch to GPUs with right computer:
        if (Input.sys.host_idn == 6) || Input.sys.host_idn == 2
            Input.inversion_tool='fft_gmres_gpu';
        end

        if Input.sys.host_idn > 0
            Input.rootfolder=Input.sys.rootfolder_list{Input.sys.host_idn};
        end        
        if ~isfield(Input,'rootfolder')
            DDA_simulation_error(error_number);
        end        
                                         
        %% particle file, DDA engine etc. related
        Input.dipole_spacing=2e-9;        % dipole spacing, now either 1e-9 or 2e-9
                      
        %% Array related parameters:
        Input.outputs = 'C_P_and_E';      % handle to control if P is also output 'C_P_and_E' , 'C_and_P' or 'C'
        %Input.outputs = 'C';      % handle to control if P is also output 'C_P_and_E' , 'C_and_P' or 'C'
        Input.save_results = true;      % handle to control if results are saved or not

        Input.sample.n=1.51;            % 1.51 refractive index of surrounding medium
        Input.n=Input.sample.n;
        
        Input.arrayN=200;                         % array size is NxN    215
        Input.lambda_peak_vec=1e-9*[515 515 515]; % plasmon resonances lambda_x, lambda_y, lambda_z
        Input.particle.tau=0.5e-15;               % used to be 0.3e-15 for silver, check how could change things 0.05e-15
        Input.particle.scatter_strength=0.04e-5; % A0=0.3e-5 1Dz: 0.05e-5!!! THIS IS HUGELY IMPORTANT !!

        % Gmres input parameters:
        Input.gmres.restart=10;
        Input.gmres.maxit=1000;
        Input.gmres.tol=2e-6;

        Input.dimensionality='3D';            % 1Dz, 3D, 1Dx
%         Input.array_type='stl_file';      % lattice type :    'square' 'hexagonal'
%         Input.sample.type='stl_file';     % 'lorentzian Au', 'ring', 'oligomer'
        
        Input.filename='test2.stl';  % Spiral_2loops2.stl
        Input.stl_filename='';  % Spiral_2loops2.stl
        
        %% Excitation field related parameters:
        % 'x-plane' = single plane wave from angle theta, i.e. kx != 0, ky=0
        % '2 x-planes' = 2 plane waves from angles +-theta, i.e. (kx,0), (-kx,0)
        % '3 x-planes' = 3 plane waves from: (kx,0) , and it's 120 deg rotated counterparts
        % '4 x-planes' = 4 plane waves from: (kx,0) , (-kx,0) , (0,ky=kx), (0,-kx)
        Input.excitationtype='x-plane';

        Input.theta=0;      % Incident angle (deg), 0 = normal incidence
        Input.phi=0;        % angle in xy-plane, 0 = along x, rot towards y-axis
        Input.psi=0;        % rotation of x'y'-plane, 0: x=x', rot towards y'-axis
        % e.g. angles (0,0,0) give y-polarized plane along z-axis
        % e.g. angles (0,0,90) give x-polarized plane along z-axis
        
        Input.N=Input.arrayN.*Input.arrayN; % number of particles in array       
    end

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

end

