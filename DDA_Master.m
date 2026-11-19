function DDA_Master(Input)
  %close all
  
    Make_predefined_Inputs()
    if isfield(Input, 'run_cancelled')
        disp('Run cancelled in the parameter GUI.')
        return
    end
    Save_function_workspace_when_debugging()
    %Make_particle_inputs()
  
    Result=struct; 
   % disp(['Scatter cross-section: ' num2str(Input.particle.scatter_strength)])
         
    Make_savename();
         

    [Result.C,Result.P,~,~] = DDA_calc_looped(Input);
      
    toc  
    % if ~isfield(Input,'NLO_modality')
        % handle to control if results are saved or not
        if  Input.save_results && ~isempty(Result.P)   % (empty for submitted jobs)
            [Result] = Save_DDA_results(Input,Result);
        end
    % end
    %% plot 2D if wanted
    % Plot_2D_period_vs_wavelength(Result(1).filename,plotting);
    % Plot_2D_angle_vs_wavelength(filename)
    
    % Function to make the savename
        function Make_savename()
            if isfield(Input, 'NLO_modality') % length(Input.NLO_modality) > 1
                 %[~,~,Input.filename] = Make_DDA_save_folder_and_name(Input,Input.NLO_modality); 
                 [~,~,Input.filename] = Make_DDA_save_folder_and_name_TS(Input,Input.NLO_modality); 
                %disp(['Filename: ' Input.filename])
            else
                %              Input.filename=['TEST'];
                if isfield(Input.particle,'width')
                    w_str=['afit_L' num2str(Input.particle.L*1e9) ...
                        '_W' num2str(Input.particle.width*1e9)];                
                else
                    w_str=['lspr' num2str(Input.particle.lambda_peak_vec(1)*1e9)];
                end
                
                % form a string variable from periodicities:
                period_str = ['P' char(strrep(join(string(round(Input.p*1e9,1)),'_'),'.','p'))];
                 
                % form string from array size parameters:
                array_str = ['N' char(join(string(Input.arrayN),'_'))];
                 
                Input.filename=['L_' Input.particle.particle_type '_' ...
                    period_str ...
                    '_Ei' num2str(Input.theta) ...
                    '_' num2str(Input.phi) ...
                    '_' num2str(Input.psi) ...
                    '_wl' num2str(round(nm2ev(max(Input.range.eV))*1e9)) ...
                    '_'   num2str(round(nm2ev(min(Input.range.eV))*1e9)) '_' ...
                    array_str '_' w_str];
                
                %disp(['Filename: ' Input.filename])
            end
        end
    
    % Function to save the function workspace when debugging
        function Save_function_workspace_when_debugging()
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

%% plotting inputs, nested function
    
    function Make_predefined_Inputs 
       %close all
        %% Parameter GUI
        Input.open_parameter_GUI = true;   % true: edit the options below in a window before the run (DDA_parameter_GUI)
        Input.use_saved_defaults = true;   % true: apply DDA_saved_defaults.m (written by the GUI's "Save defaults")
                                           % after the values below; delete that file to use only these values
        Input.favorite = "";               % name of a saved favorite (Favorites tab of the GUI) to apply
                                           % after the saved defaults; "" = none

        %% array inputs
        Input.save_message = "";%"Outerwall y 30nm (1)"; 
        use_eV = 1;
        Input.n = 1.4665; % not in use if "substrate" is used
        Input.p = [700,700,100]*1e-9;  %Input.p = [435,690,500]*1e-9;  %435nm: period of the measured DFB pixels
        Input.arrayN=[20 20 1]; %80µm pixel has 184 period if 435nm is the period
        Input.dimensionality="3D";            % Dipolemoment dimensionality, '1Dx', '3D', '2D'etc for point like dipoles,
                                              % "cyl_TM" infinite cylinders along y, field along the cylinders (scalar, N unknowns),
                                              % "cyl_TE" the same cylinders, field in the xz plane (p_x and p_z, 2N unknowns).
                                              % Cylinder arrays run along x: arrayN = [Nx 1 1].

        %% Layered background for the cylinder models (cyl_TM, cyl_TE), see dda_layers
        % The cylinders sit in the "host layer", whose index is the host index above
        % (Input.n with Input.substrate = 'const'). Any number of layers below and
        % above it, listed from the host layer outwards; the last index of each list
        % is the half-space, so it has no thickness.
        % Example, air / 1.18 um film / 2 um MY133 (1.33) / glass, grating ridges as
        % cylinders of index 1.33 (Input.Mie_glass_n) lying on the lower interface:
        %   Input.n = 1.4665; Input.substrate = 'const'; Input.sample.type = "glass_Mie";
        %   Input.layers.host_thickness = 1.18e-6;  Input.layers.z_cyl = Input.particle.r;
        %   Input.layers.below_n = [1.33 1.52];     Input.layers.below_d = 2e-6;
        %   Input.layers.above_n = 1;               Input.layers.above_d = [];
        Input.layers.enable         = false;       % true: use the layered background
        Input.layers.host_thickness = 1.18e-6;     % thickness of the layer with the cylinders [m]
        Input.layers.z_cyl          = 100e-9;      % height of the cylinder axes above its lower interface [m]
        Input.layers.below_n        = [1.33 1.52]; % indices below, nearest first; last = substrate half-space
        Input.layers.below_d        = 2e-6;        % thicknesses below [m], one less than below_n
        Input.layers.above_n        = 1;           % indices above, nearest first; last = cover half-space
        Input.layers.above_d        = [];          % thicknesses above [m], one less than above_n
        Input.particle.cyl_alpha_legacy = false;   % true: the cyl_TM polarizability used before the fix
                                                   % (too small; only to reproduce old results)


        %% Boundaries of the interaction matrix
        Input.mirror_A = [0,0,0,0]; % mirrors at [x=0, x=max(x), y=0, y=max(y)]: 0 none, +1 or -1.
                                    % The planes go through the outermost particles: those are their own
                                    % images (a single column between x mirrors is not changed by them).
                                    % 2D/3D: +1 field mirror-symmetric (tangential p even, normal p odd; PMC-like),
                                    %        -1 antisymmetric (tangential odd, normal even; PEC-like).
                                    % 1Dx/1Dy/1Dz: p(-x) = +-p(x) for the one component.
        Input.periodic_boundary = [0,0]; %[x,y], Must not overlap with .mirror_A
        % Input.periodic_cell   = [700e-9 0];  % cell size [m]; needed for a single row/column
        % Input.periodic_Ncells = [10 10];     % image cells on EACH side (0 = none, 10 = 21 cells)
        % -1 for p(-x) = -p(x). Boundaries: [x=0, x=L, y=0, y=L]

        %% Wavelengths and range
        wl = 1005.5:.5:1007;
        wl = linspace(1042,1068,100);
        Input.range.imag_parts = linspace(-30, -0, 10)*1e-3; %eV, set to 0 for real analysis
        %Input.range.imag_parts = linspace(0,0,1)*1e-3; %eV, set to 0 for real analysis
        Input.range.real_parts = linspace(1.1,1.5,40);
       % Input.range.real_parts = linspace(1240/1100,1240/1030,50);
        

        %Input.range.eV = wl_to_eV./Input.range.wl;
       % Input.range.eV=reshape(nm2ev(Input.range.wl)+Input.range.maginary_parts.'*1i,1,length(Input.range.wl)*length(Input.range.imaginary_parts));
       wl_to_eV=1239.84187e-9;
       if use_eV
            Input.range.eV =  reshape(Input.range.real_parts +Input.range.imag_parts.'*1i,1,length(Input.range.real_parts)*length(Input.range.imag_parts));
            Input.range.wl = wl_to_eV./Input.range.eV;
       else
           Input.range.wl = wl*1e-9;
            
       end

        c = 2.99792458e8;
        %Input.range.k = sqrt((eV2rads(Input.range.eV)).^2./c.^2);
        Input.range.k = 2*pi./Input.range.wl;
        Input.range.k_N = length(Input.range.k);      
          %% Mie resonance parameters
        Input.particle.r = 65e-9;   %"limit" for dipole approximation
        %Input.particle.r = 110e-9; %Glass cylinders 


        Input.particle.anisotropic = false;
        Input.particle.rx = 65e-9;
        Input.particle.ry = 100e-9;
        Input.particle.rz = 10e-9;
        %% Particle types and arguments, substrate
        Input.substrate = 'BK7';
        %Input.substrate = 'const';
        % Particle model:
        %   Mie theory (Input.particle.r):  "Au_Drude_Mie" | "Au_Mie" | "glass_Mie"
        %   MJH Lorentzian dipoles:         "Lorentzian"  (see "MJH particle arguments" below)
        Input.sample.type="Au_Drude_Mie"; %"glass_Mie", "Au_Drude_Mie", "Lorentzian"
        Input.Mie_glass_n=1.476;Input.particle_gain=0.;
        
        Input.array_type = '1Dsquare_halfshifts' ;      % lattice type : 'PC_mirrors' 'double_rectangular' '1Dsquare_halfshifts'   'square' 'hexagonal' 'double_retangular'

        %Input.array_type = 'PC_mirror_reflection_test';
        %Input.array_type = 'PC_mirrors';
        Input.array_shifts_y = 0;Input.array_shifts_x = 0; %also used in outerwall, HG_wall and Inverted_gaussian!!!!!!
        Input.PC_mirror_thickness = 0; %NPs at each side. Total number of particles still determined by Input.arrayN. PC_mirrors dont yet work with array_shifts
        Input.PC_mirror_p = [700e-9,350e-9];
        %Input.PC_mirror_p = [700e-9, 700e-9/4];
        Input.PC_mirror_offset = [0e-9, 0e-9];

         if Input.array_type == "double_rectangular"
             Input.nmbofmodes = Input.arrayN(1)*Input.arrayN(2)*2;
         else
             Input.nmbofmodes = Input.arrayN(1)*Input.arrayN(2);
         end
        Input.second_array_offset = [Input.p(1)/4, Input.p(2)/4,0];

       % Input.array_type = 'Outerwall'; 
        Input.Outerwallxy = 2; % for wall on: both x and y: 0, only x: 1, only y: 2
        %Input.array_type = "HG_wall";
        %Input.array_type = "Gaussian"; %inverted 1D gaussian of particle sizes in y direction. r+"wall_increase_in_radius" at edges and r at center
        Input.HG_wall_mode = [0,0];
        Input.wall_increase_in_radius=0e-9;
        Input.wall_thickness = 1; %nmb of particles
        Input.Gaussian_sigma = Input.arrayN(1)*Input.p(1)/4*1.5;
       
        %Input.array_type='shifted_diagonals';  %not worthit
        %% Excitation field related parameters:

        % 'x-plane' = single plane wave from angle theta, i.e. kx != 0, ky=0
        % '2 x-planes' = 2 plane waves from angles +-theta, i.e. (kx,0), (-kx,0)
        % '3 x-planes' = 3 plane waves from: (kx,0) , and it's 120 deg rotated counterparts
        % '4 x-planes' = 4 plane waves from: (kx,0) , (-kx,0) , (0,ky=kx), (0,-kx)
        Input.polarization ='linear'; %'RC, 'LC', 'linear'.. right and left circular only for normal incidence

        Input.theta=0.;      % Incident angle (deg), 0 = normal incidence
        Input.phi=0;        % angle in xy-plane, 0 = along x, rot towards y-axis.
        % Phi rotates also polarization. 90 deg difference between phi and psi gives x polarized. If using 1D simulations, remember to change both
        Input.psi= 90; %30 gives best gain-assisted Q-factor (for slr in p_y)?
        % e.g. angles (0,0,0) give y-polarized plane along z-axis
        % e.g. angles (0,0,90) give x-polarized plane along z-axis
         % e.g. angles (0,90,0) give x-polarized plane along z-axis
        %Great Q factor:
        %Gain 0.00525, psi = 30, arrayN =65,65; p = 500 1060, radius = 70nm


        % Field profile: "none" (plane wave) or "Gaussian": a Hermite-Gaussian beam
        % HG(Hermite_order_x, Hermite_order_y) centred on the array, waist w0 at a
        % distance Gaussian_z from the array (z > 0: waist beyond the array).
        % Cylinders: 1D beam in x (Hermite_order_y not used).
        Input.preset_field = "none";
        Input.Hermite_order_x = 0;  Input.Hermite_order_y = 0;
        FWHM = 25e-6;                              % intensity FWHM at the waist
        Input.Gaussian_w0 = FWHM/sqrt(2*log(2));   % waist radius (1/e^2 intensity)
        Input.Gaussian_z  = 0;
        %% Other parameters

        Input.gain = 0.0e-3;
        
        %not important, obsolete tests
        Input.Hermitian = false;
       
        %Input.Gaussian.theta = asin(3/f);

        Input.try2Das1D = 0;
        Input.nearest_neighbour =false; %only nearest interactions taken into account(in fft_gmres)

        Input.diag_shift = 0e-9; %if Input.array_type='shifted_diagonals', dictates how much every other particle is shifted diagonally 
       %% MJH particle arguments (used when Input.sample.type = "Lorentzian")
        % Per axis c: alpha_c(w) = A0_c/(w0_c - w - i*gamma_c), w0 = 2*pi*c/lambda_peak,
        % gamma = 1/(2*pi*tau). Converted to the solver's SI units
        % (alpha_SI = 4*pi*eps0*n_host^2*alpha) by Lorentzian_polarizability_DDA.
        % Each of the three can be a scalar, [x y z], or N x 3 (one row per particle,
        % e.g. for particles that are not identical).
        Input.particle.scatter_strength = 0.01e-5;          % A0; 0.3e-5, 1Dz: 0.05e-5. THIS IS HUGELY IMPORTANT
        %Input.particle.scatter_strength = [1 1 1]*1e-5/(4*pi);
        Input.particle.tau = [1 4 4]*1e-15;                 % lifetimes [s] (x,y,z); was [2 0.2 2]*1e-15
        Input.particle.lambda_peak_vec = 1e-9*[780 330 1000];  % resonance wavelengths (x,y,z)
        Input.particle.lorentz_units = "gaussian";          % "gaussian": A0 values as in the old MJH code
                                                            % "SI": A0 already gives alpha in C m^2/V
        Input.particle.lorentz_correction = "none";         % "none" (as before) | "radiative" (energy
                                                            % conserving) | "MLWA" (+ depolarization)
        Input.particle.lorentz_mlwa_a = 100e-9;             % particle size used by "MLWA"
        Input.particle.particle_type = 'Au_L';              % only used in Save_DDA_results file names
        %% General arguments, other MJH plotting arguments, 
        warning('off')
        Input.Eloc=true;
        Input.debug=false;
        Input.project='DDA';
        Input.outputs = 'C_P_and_E';      % handle to control if P is also output 'C_P_and_E' , 'C_and_P' or 'C'
        %Input.outputs = 'C';      % handle to control if P is also output 'C_P_and_E' , 'C_and_P' or 'C'
        Input.save_results = false;      % handle to control if results are saved or not
        
        
        %% =====================================================================
        %% Computing options
        %% =====================================================================
        %% Solver
        Input.inversion_tool='eig';      % linsolve, mldivide , gmres_gpu, 'fft_gmres', 'gmres'
        %Input.inversion_tool = 'fft_gmres';
        %Input.inversion_tool = 'gmres';
        %Input.inversion_tool='inv';
        Input.useparallel = 1;           % read only by the old solvers; see Input.parallel below

        % Gmres input parameters:
        Input.gmres.restart=100;Input.gmres.maxit=100;Input.gmres.tol=2e-6;

        %% Where the solve runs
        % "none"       : in this MATLAB, serial; MATLAB is busy until done
        % "threads"    : in this MATLAB, local thread pool (Input.threads.*); busy until done
        % "background" : local pool via parfeval, MATLAB stays free (Input.background.*).
        %                Each run uses one worker; several runs go in parallel.
        % "cluster"    : batch job via MATLAB Parallel Server (Input.cluster.*), MATLAB
        %                stays free unless Input.cluster.wait. Input.cluster.name =
        %                "Processes" runs the same kind of job on this computer, with
        %                its own pool (one big run in the background).
        % Submitted runs: collect_dda_jobs fetches and saves the finished ones,
        % dda_jobs lists, cancels or removes them.
        Input.parallel = "threads";
        Input.threads.num_workers = [];        % [] = all cores
        Input.threads.keep_pool   = true;      % keep the pool open between runs

        Input.background.pool        = "Processes";  % "Processes" (robust) or "Threads"
        Input.background.num_workers = [];           % [] = MATLAB's default pool size

        % Cluster job. Cores reserved = (num_workers + 1) * num_threads,
        % because batch uses one extra worker to run the client code.
        Input.cluster.name        = "Roihu R2026b";  % parcluster profile ("Processes" = this computer)
        Input.cluster.wait        = false;           % true: block until done, then save and plot
        Input.cluster.wall_time   = "0:03:00";       % h:mm:ss, job is killed after this. 3Dx20x20 solve 2.9s / single thread worker
        Input.cluster.num_threads = 1;               % threads per worker (BLAS/eig)
        Input.cluster.num_workers = 10;              % pool workers running the parfor
        Input.cluster.partition   = "";              % "" = profile default
        Input.cluster.project     = "";              % CSC project (Slurm AccountName), e.g. "project_2001234"; required on Roihu.
                                                     % Set yours in the GUI and press "Save defaults" (DDA_saved_defaults.m, git-ignored)
        Input.cluster.mem_per_cpu = "";              % e.g. "4g"; "" = profile default
        Input.cluster.email       = "";              % Slurm notifications
        Input.cluster.submit_args = "";              % extra sbatch arguments
        Input.cluster.additional_properties = struct(); % any other AdditionalProperties entry
        % If your plugin names a property differently, map it here, e.g.
        % Input.cluster.property_names.partition = 'QueueName';
        Input.cluster.property_names = struct();
        Input.cluster.current_folder       = ".";
        Input.cluster.auto_add_client_path = false;
        Input.cluster.auto_attach_files    = true;   % send the needed .m files with the job
        Input.cluster.attached_files       = {};
        Input.cluster.additional_paths     = {};     % code folders that already exist on the cluster
        Input.cluster.delete_job   = true;          % delete after fetching (failed jobs are kept)
        Input.cluster.save_profile = false;         % store these settings in the profile

        Input.jobs.plot_on_collect = true;    % plot submitted runs when collect_dda_jobs saves them

        %% Eigenmode tracking (Track_eigenmodes, called by Calculate_DDA_eigenmodes)
        % Gives every eigenvector a (mx,my,pol) label and follows it over frequency.
        % Needed for the eigenmode plots (C.extsall, C.coefs_ordered, C.P_eigs).
        Input.track.enable         = "auto";    % "auto" = real-frequency runs only | true | false
        Input.track.method         = "hybrid";  % "basis" | "overlap" | "hybrid"
        Input.track.overlap_weight = 0.5;       % weight of eigenvector overlap in "hybrid"
        Input.track.basis          = "auto";    % "auto" (Fourier along periodic directions) | "sine" | "fourier"
                                                % | "k" (plane wave that couples best, as the complex plots' "k" labels)
        Input.track.k_step         = false;     % "k": undo the pi steps of half-shifted lattices first
        Input.track.store_P_eigs   = true;      % per-mode dipole moments, Ndof x Nk each
        Input.track.min_confidence = 0.5;       % basis score below this counts as unclear
        Input.track.max_memory_GB  = 8;         % guard: all eigenvectors are kept while tracking
        Input.track.verbose        = true;

        Input.degeneracy_tol = 1e29;
        Input.P_method = "direct";               % driven response in the eigen solver: "direct" P = A\E (exact),
                                                 % "eigenmodes" P = V*(c./lambda), the sum over eigenmodes            % |lambda - lambda_min| tolerance for C.degeneracy
        Input.use_legacy_eig_solvers = false;   % true = old per-dimension files (for validation);
                                                % their results are saved with ' legacy' in the name
        Input.recompute = true;                % true = ignore an existing result file and overwrite it
        %% =====================================================================
        %% Saving and plotting
        %% =====================================================================
       %% Saving real-frequency results (save_real_freq_results)
       Input.save.P_eigs      = true;   % per-mode dipole moments, needed for the mode plots
       Input.save.max_file_GB = 4;      % above this C.P_eigs is dropped (with a warning)
       Input.save.complex_all = false;  % complex frequencies: false keeps only what the plots use (smallest
                                        % eigenvalue, its eigenvector, degeneracy); true keeps all eigenvalues,
                                        % coefficients, P and the incident field (about 8x larger)
       %% Plotting (real_freqs_plots_from_file, complex_freqs_plots_from_file)
       % Results are always saved; these only decide what is plotted. To re-plot a
       % saved run with other options: real_freqs_plots_from_file(file, plotting)
       % Wavelengths in nm, eigenmode labels 0-based [mx my] or [mx my pol].
       Input.plotting.enable        = true;      % plot after the run
       Input.plotting.export_folder = "";        % "" = no export, else PNGs of every figure go here
       Input.plotting.docked        = true;      % true: all figures as tabs in one window; false: separate windows

       % Transmission, |p| vs position and wavelength, profiles and maps
       Input.plotting.dipoles.enable      = true;
       Input.plotting.dipoles.modes       = [0 0; 0 1]; % eigenmodes shown on their own (eig)
       Input.plotting.dipoles.n_strongest = 0;       % + the N modes most excited at the extinction peak
       Input.plotting.dipoles.profile_wls = 1060.4;  % profiles through the centre column (+ the peak)
       Input.plotting.dipoles.map_wl      = [];      % |p| map wavelength, [] = first wavelength
       Input.plotting.dipoles.surf_wls    = [];      % |p| surfaces coloured by phase, e.g. [1053 1066]

       Input.plotting.Q.enable = true;              % Q-factor fit of the extinction
       Input.plotting.Q.wl     = [];                % [] = extinction peak
       Input.plotting.Q.guess  = 30;

       Input.plotting.incidence.enable = true;      % incident field and its angular spectrum
       Input.plotting.incidence.wl     = 800;

       Input.plotting.apertures.enable = false;     % extinction through central apertures
       Input.plotting.apertures.sizes  = [1 10 20 30 40 60 80];   % particles per side

       Input.plotting.emitted.enable = false;
       Input.plotting.emitted.types  = "fraunhofer";  % any of "fraunhofer" "transmitted" "far_field" "quiver_ball"
       Input.plotting.emitted.wls    = 630;           % [] = extinction peak
       Input.plotting.emitted.angle_lim = 2;          % deg, axes of fraunhofer/transmitted
       Input.plotting.emitted.fraunhofer_quantity = "P/alpha";  % "P/alpha" (as before) or "P"
       Input.plotting.emitted.anglex = linspace(-0.03,0.03,20); % rad, far_field view
       Input.plotting.emitted.angley = [];            % [] = same as anglex

       Input.plotting.coefficients.enable  = true;  % weight of the most excited eigenmodes (eig)
       Input.plotting.coefficients.n_modes = 10;
       Input.plotting.eig_band.enable      = true;  % (mx,0) and (0,my) mode extinctions vs energy (eig)

       % Complex-frequency plots (complex_freqs_plots_from_file); saved with the
       % results, so a file can be re-plotted with other options:
       % complex_freqs_plots_from_file(file, struct('quiver', false))
       Input.plotting.complex.full_range       = true;        % axes over the whole computed range
       Input.plotting.complex.real_axis        = "eV";        % real axis: "eV" (Re of the frequency) or "nm" (wavelength)
       Input.plotting.complex.re_lims          = [1.12 1.18]; % eV,  used when full_range = false
       Input.plotting.complex.im_lims          = [25 30];     % meV (-Im), used when full_range = false
       Input.plotting.complex.modes            = true;        % plot the modes at the minima of |lambda|
       Input.plotting.complex.tiles_per_figure = 20;
       Input.plotting.complex.max_fft_order    = 1000;        % skip modes with a higher spatial order
       Input.plotting.complex.quiver           = true;        % 2D/3D modes as arrows (needs slanCM)
       Input.plotting.complex.imagesc_is_z     = true;        % 3D: colour map shows p_z
       Input.plotting.complex.separate_pols    = false;       % 2D/3D: one tile per polarization
       Input.plotting.complex.surf_1D          = true;        % 1D-dimensional modes as surfaces
       Input.plotting.complex.labels           = "sine";      % mode labels: "sine" (basis of the mode tracking),
                                                              % "k" (plane wave that couples best), "nodes" (sign changes),
                                                              % "fft" (as before)
       Input.plotting.complex.k_step           = false;       % "k": undo the pi steps of half-shifted lattices first
       Input.plotting.complex.basis            = "auto";      % "auto" (Fourier along periodic directions) | "sine" | "fourier"
       Input.plotting.complex.min_confidence   = 0.25;        % sine/nodes labels below this get a '?'
       Input.plotting.complex.label_map        = false;       % map of the nearest mode's number over the complex plane
       Input.plotting.complex.plot_fft         = true;        % show the spatial FFT used for labels ("fft" labels only)
       Input.plotting.complex.nufft            = false;       % non-uniform FFT for the labels ("fft" labels only)
        %% Saved defaults: DDA_saved_defaults.m overrides the values above
        code_defaults = Input;
        if Input.use_saved_defaults && exist('DDA_saved_defaults', 'file') == 2
            Input = DDA_saved_defaults(Input);
            Input = dda_update_derived_inputs(Input, code_defaults, use_eV, wl_to_eV);
        end

        %% Favorite: a saved input profile (Input.favorite), e.g. for scripted runs
        if isfield(Input, 'favorite') && strlength(string(Input.favorite)) > 0
            before_fav = Input;
            Input = dda_favorites('apply', Input, Input.favorite);
            Input = dda_update_derived_inputs(Input, before_fav, use_eV, wl_to_eV);
        end

        %% Parameter GUI: edit the options above (Input.open_parameter_GUI)
        if Input.open_parameter_GUI
            before = Input;
            [Input, gui_ok] = DDA_parameter_GUI(Input, code_defaults);
            if ~gui_ok
                Input.run_cancelled = true;
                return
            end
            Input = dda_update_derived_inputs(Input, before, use_eV, wl_to_eV);
        end

        %% system arguments
        if Input.gain ~= 0

            disp("It is recommended to use complex frequencies with Input.range.imaginary_parts as opposed to Input.gain ")
            
        end
        % Make rootfolder for the calculations by taking the current folder as the master:
        p = mfilename('fullpath');
        [pathstr,name,~]=fileparts(p);
        Input.code_rootfolder=pathstr;
        Input.m_name=name;
        Input.stl_folder=[Input.code_rootfolder 'stl_samples\'];
            
        % Host-specific folders (optional): Get_cpu_list_and_folder_hierarchy_TS.m
        % holds machine names and local paths, so it is kept out of the
        % repository (.gitignore). Without it the DDA folder itself is the root.
        if exist('Get_cpu_list_and_folder_hierarchy_TS', 'file') == 2
            [Input.sys] = Get_cpu_list_and_folder_hierarchy_TS();
        else
            Input.sys = struct('host_idn', 0);
        end
        Input.date=date;
        
        % switch to GPUs with right computer:
        if (Input.sys.host_idn == 6) || Input.sys.host_idn == 2
            Input.inversion_tool='fft_gmres_gpu';
        end

        if Input.sys.host_idn > 0
            Input.rootfolder=Input.sys.rootfolder_list{Input.sys.host_idn};
        else
            Input.rootfolder=[Input.code_rootfolder filesep];
        end        
        %array from file arguments
         Input.filename='test2.stl';  % Spiral_2loops2.stl
        Input.stl_filename='';  % Spiral_2loops2.stl 
%         Input.array_type='stl_file';      % lattice type :    'square' 'hexagonal'
%         Input.sample.type='stl_file';     % 'lorentzian Au', 'ring', 'oligomer'
        %% particle file, DDA engine etc. related, OBSOLITE

        %Input.dipole_spacing=2e-9;        % dipole spacing, now either 1e-9 or 2e-9
                      
    end

    
end

function [C_out,P_out,Field_out,Array_out,Input] = DDA_calc_looped(Input)
    %% Solve polarization and cross-secttions for arrays, using DDA. Can use 
    % either 1D case, where only N x N matrices need to be inverted, or full 3D
    % case where 3N x 3N matrices are operated on.
    %  
    tic 
    Input.toc1=toc;
    
    %% 1. Additional parameters calculated from inputs:
    if length(Input.arrayN) >2
        Input.N=Input.arrayN(1).*Input.arrayN(2)*Input.arrayN(3); % number of particles in array
    else
        Input.N=Input.arrayN.*Input.arrayN; % number of particles in array
    end
    
    %% 2. Make Array struct containing info of the scatterers, locations etc.
    [Array,Input]=Make_generic_Array(Input);
    [Array,Input]=dda_apply_array_edit(Array,Input);   % layout edited in the array editor (GUI), if any
    mirror_notice(Input, Array);
    
    %% 3. Make single particle polarizabilities, linear and NLO:
    [Array]=NLO_DDA_Make_3D_polarizability(Array,Input);
    %% 4. Solve, save and plot
    % Both pipelines: name the result file from the inputs; if the file exists
    % (and holds valid results) plot it, otherwise solve, save, then plot it.
    is_complex = any(Input.range.imag_parts);
    if Input.array_type == "PC_mirror_reflection_test"
        PC_mirror_reflection_test(Input,Array);
        C_out = []; P_out = []; Field_out = []; Array_out = Array;
        return
    end
    if is_complex && Input.inversion_tool ~= "eig"
        error('Complex frequencies need Input.inversion_tool = ''eig''.');
    end

    if is_complex
        kind = "complex";
    else
        kind = "real";
    end
    filename = construct_filename(Input, kind);
    if Input.use_legacy_eig_solvers && Input.inversion_tool == "eig"
        % construct_filename does not record the solver; keep both results apart
        filename = replace(filename, ".mat", " legacy.mat");
    end
    fullpath = fullfile(dda_results_folder(kind), filename);

    info = struct('fullpath', fullpath, 'filename', filename, 'kind', kind, 'Input', Input);
    mode = "none";
    if isfield(Input,'parallel'), mode = string(Input.parallel); end
    legacy = Input.use_legacy_eig_solvers && Input.inversion_tool == "eig";
    C_out = []; P_out = []; Field_out = []; Array_out = Array;

    % (a) A valid result file exists: plot it.
    if isfile(fullpath) && ~Input.recompute && ~cached_result_is_invalid(fullpath)
        fprintf("Already calculated, plotting from file (set Input.recompute = true to redo):\n  %s\n", fullpath)
        S = load(fullpath, "C");  C_out = S.C;
        if Input.plotting.enable
            if is_complex, complex_freqs_plots_from_file(filename)
            else,          real_freqs_plots_from_file(fullpath)
            end
        end
        return
    end

    % (b) This run is already submitted and has not failed: do not submit it twice.
    pending = dda_job_registry("find", fullpath);
    if ~isempty(pending)
        state = dda_job_registry("state", pending);
        if any(state == ["queued","running","finished"])
            fprintf("Already submitted (%s); collect_dda_jobs fetches it when it is done:\n  %s\n", state, filename)
            C_out = struct('job', pending.key, 'state', state);
            return
        end
        dda_job_registry("remove", pending);          % failed or lost: submit again
    end

    % (c) Solve.
    switch mode
        case {"cluster","background"}                 % submitted; MATLAB stays free
            if legacy
                error('The legacy eigen solvers (validation) run only with Input.parallel = "none" or "threads".');
            end
            entry = submit_dda_job(Input, Array, fullpath, filename, kind);
            if mode == "cluster" && isfield(Input.cluster,'wait') && Input.cluster.wait
                dda_job_registry("wait", entry);
                if dda_job_registry("state", entry) ~= "finished"
                    [msg, tail] = dda_job_registry("message", entry);
                    error('Cluster job %d failed: %s\n%s\nThe job is kept for inspection; see dda_jobs.', ...
                          entry.job_id, msg, tail);
                end
                R = dda_job_registry("fetch", entry);
                finalize_dda_result(R, info, Input.plotting.enable);
                dda_job_registry("remove", entry);
            else
                fprintf("Run collect_dda_jobs to fetch finished jobs, dda_jobs to list them.\n")
                C_out = struct('job', entry.key, 'state', "submitted");
                return
            end
        otherwise                                     % "none" / "threads": in this MATLAB
            if legacy
                [R.P, R.C, R.Field, R.Array] = legacy_eig_dispatch(Input, Array);
            else
                R = dda_solve_job(Input, Array, mode);
            end
            finalize_dda_result(R, info, Input.plotting.enable);
    end
    C_out = R.C;  P_out = R.P;  Field_out = R.Field;  Array_out = R.Array;

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

function bad = cached_result_is_invalid(filename)
    % Treat an existing result file as missing when
    %  - it holds only NaN (the old 3D complex solver in cluster mode never
    %    fetched its batch results), or
    %  - it is a 2D/3D eigen-solver result from before the self-term fix in
    %    Make_int_matrix_A_SI (no C.solver_version >= 2): those matrices used
    %    an element-wise 1./alpha and are wrong.
    bad = false;  why = "";
    try
        S = load(filename, "C", "Input");
        if isfield(S.C,'min_eig_values')
            bad = all(isnan(S.C.min_eig_values(:)));
        else
            bad = ~isfield(S.C,'ext') || all(isnan(S.C.ext(:)));
        end
        if bad, why = "no valid results"; end
        if ~bad && isfield(S,'Input') && string(S.Input.dimensionality) == "cyl_TM" ...
                && string(S.Input.inversion_tool) == "eig" ...
                && ~(isfield(S.C,'solver_version') && S.C.solver_version >= 3) ...
                && ~(isfield(S.Input.particle,'cyl_alpha_legacy') && S.Input.particle.cyl_alpha_legacy)
            bad = true;  why = "cyl_TM result from before the cylinder polarizability fix";
        end
        if ~bad && isfield(S,'Input') && isfield(S.Input,'mirror_A') && any(S.Input.mirror_A) ...
                && ~(isfield(S.C,'solver_version') && S.C.solver_version >= 4)
            bad = true;  why = "mirror result from before the mirror fix (particles on the plane counted twice)";
        end
        if ~bad && isfield(S,'Input') && any(string(S.Input.dimensionality) == ["2D","3D"]) ...
                && string(S.Input.inversion_tool) == "eig" ...
                && ~(isfield(S.C,'solver_version') && S.C.solver_version >= 2)
            bad = true;  why = "2D/3D result from before the self-term fix";
        end
    catch
        bad = true;  why = "file could not be read";
    end
    if bad
        fprintf("Existing file ignored (%s), recalculating: %s\n", why, filename);
    end
end

function mirror_notice(Input, Array)
    % The mirror planes go through the outermost particles (Make_int_matrix_A_SI).
    % Say so when a mirror has nothing to reflect, e.g. x mirrors on a 1 x N column.
    if ~isfield(Input,'mirror_A') || ~any(Input.mirror_A), return; end
    r = double(Array.r_i);
    tol = 1e-3*min(Input.p(1:2));
    m = zeros(1,4); m(1:numel(Input.mirror_A)) = Input.mirror_A;
    planes = [0, max(r(:,1)), 0, max(r(:,2))];
    names  = ["x = 0", "x = max(x)", "y = 0", "y = max(y)"];
    for j = find(m)
        d = 1 + (j > 2);
        if all(abs(r(:,d) - planes(j)) < tol)
            if m(j) > 0
                fprintf(['Mirror at %s: every particle lies on this plane, so it changes nothing ' ...
                         '(the particles are their own images).\n'], names(j));
            else
                fprintf(['Mirror at %s (-1): every particle lies on this plane, so the symmetry ' ...
                         'forbids the components that are odd there; they are decoupled.\n'], names(j));
            end
        end
    end
end

function [P,C,Field,Array_out] = legacy_eig_dispatch(Input,Array)
    % The old per-dimension solvers, kept only to validate Calculate_DDA_eigenmodes.
    % They read flags that are no longer set in Make_predefined_Inputs; give them
    % their old values here (no plots from inside the solvers).
    old = struct('plotall',0, 'make_eig_video',0, 'plot_coefficients',0, ...
                 'plot_ffts',0, 'plot_modes1',0, 'plot_selected_modes',0, ...
                 'fft_mode_search', double(Input.dimensionality ~= "2D"), ...
                 'mode_overlap_search', 2*double(Input.dimensionality == "2D"));
    f = fieldnames(old);
    for n = 1:numel(f)
        if ~isfield(Input, f{n}), Input.(f{n}) = old.(f{n}); end
    end
    Array = dda_expand_alpha(Array);      % the old solvers read the full per-particle Array.alpha
    is_complex = any(Input.range.imag_parts);
    switch string(Input.dimensionality)
        case "3D"
            if is_complex, [P,C,Field,Array_out] = Calculate_3D_DDA_complex_freqs_eigenmodes(Input,Array);
            else,       [P,C,Field,Array_out] = Calculate_3D_DDA_generic_eigenmodes(Input,Array); end
        case "2D"
            if is_complex, [P,C,Field,Array_out] = Calculate_2D_DDA_complex_freqs_eigenmodes(Input,Array);
            else,       [P,C,Field,Array_out] = Calculate_2D_DDA_generic_eigenmodes(Input,Array); end
        case "cyl_TE"
            error('There is no legacy cyl_TE solver; use Calculate_DDA_eigenmodes.');
        case "cyl_TM"
            if is_complex, [P,C,Field,Array_out] = Calculate_cyl_TM_DDA_complex_freqs_eigenmodes(Input,Array);
            else, error('There is no legacy real-frequency cyl_TM solver; use Calculate_DDA_eigenmodes.'); end
        otherwise
            if is_complex, [P,C,Field,Array_out] = Calculate_1D_DDA_complex_freqs_eigenmodes(Input,Array);
            else,       [P,C,Field,Array_out] = Calculate_1D_DDA_generic_eigenmodes(Input,Array); end
    end
    C.solver_version = 2;    % the JTK builders they use now invert the self term correctly
end
