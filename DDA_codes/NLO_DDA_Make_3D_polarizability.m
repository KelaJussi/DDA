function [Array]=NLO_DDA_Make_3D_polarizability(Array,Input)
% Make particle polarizabilities a_i in to Array struct
% Load polarizabilities, which already contain the wavelengths
% outputs:      n = refractive index
%               k = wavenumber (in the medium)
%             k_N = the size of k
%           alpha = 3x3 polarizability matrix

Save_function_workspace_when_debugging()

%% 1 .particle polarizabilities folder structures
% r = root folder,  s = sample folder, w = r+s = whole folder structure 
Array.alpha_r_folder='Polarizabilities/';
Array.alpha_w_folder=[Array.alpha_r_folder Input.sample.type '/']; %  
                  
% DEFINE INPUT WAVELENGTHS ALREADY USING eV SCALE
Array.eV=Input.range.eV;
wl_to_eV=1239.84187e-9;
Array.wl=wl_to_eV./Array.eV;
Array.nm2eV=@(wl) wl_to_eV./wl;
Array.eV2nm=@(eV) wl_to_eV./eV;

%% 2. FORM SINGLE PARTICLE POLARIZABILITY
switch Input.sample.type
    case {"Lorentzian", "lorentzian", "lorentzian Au", "nonidentical"}
        % MJH Lorentzian scatterers (Input.particle.scatter_strength, tau,
        % lambda_peak_vec), converted to the solver's SI units. Different
        % particles: give the parameters per particle (N x 3).
        [Array] = Lorentzian_polarizability_DDA(Array,Input);
    case {'oligomer','Oligomer','random_surface','stiched_array','data_fit','stl_file','ldr','cm','rr'}
        % Old MJH particle sources: Make_Lorentzian_scatterer (oligomers),
        % Make_fitted_scatterer, Polarizability_from_permittivity. They give
        % non-SI polarizabilities, set neither Array.k nor Array.n, and need
        % data files that are not in the repository.
        error('NLO_DDA_Make_3D_polarizability:unsupported', ...
            ['Input.sample.type = "%s" is not supported by the SI solvers. Use "Lorentzian" ' ...
             '(per-particle parameters replace oligomers) or a Mie type.'], Input.sample.type);
    case {'Au_Mie', "glass_Mie","Au_Drude_Mie"}
        if isfield(Array, 'radius')
            % radii edited in the array editor: one Mie computation per distinct radius
            if any(string(Input.array_type) == ["Outerwall","HG_wall","Gaussian","Inverted_Gaussian"])
                error('Edited radii are not supported for %s lattices (they set their own radii).', Input.array_type);
            end
            Array = mie_per_radius(Array, Input);
        elseif any(string(Input.array_type) == ["Outerwall","HG_wall","Gaussian","Inverted_Gaussian"])
            % radii differ between particles: per particle, then one entry per distinct particle
            [Array] = Mie_polarizability_DDA(Array,Input);
            Array = dda_compact_alpha(Array);
        else
            % identical particles: compute one and share it (N times less memory)
            Array = mie_one_type(Array, Input);
        end
        %Array.alpha = eps0*Array.alpha;
    case {"Au_Drude_MLWA"}
        [Array] = MLWA_polarizability_DDA(Array,Input);
    otherwise
        error('NLO_DDA_Make_3D_polarizability:type', 'Unknown Input.sample.type "%s".', Input.sample.type);
end

% figure
% plot_complex(Array.wl,squeeze(Array.alpha(1,1,1,:)))

%% 3. put variables into GPU memory if one is in use
%{
if strncmpi('upg',rot90(Input.inversion_tool,2),3) % inversion tool uses gpu
    % variables which are used and need to be transferred to GPU mem
    Array.r_i = gpuArray(Array.r_i);
    Array.k = gpuArray(Input.sample.n*2*pi./Array.wl);
else
    Array.k = Input.sample.n*2*pi./Array.wl;    
end
Array.k_N=size(Array.k,1);
%}
%% 4. Make the nonlinear polarizabilities, beta and gamma:
% Save_function_workspace_when_debugging()

 [Array] = Make_beta_Cinfv(Array,Input);
% [Array] = Make_CTHG_beta_Cinfv(Array,Input);
if Input.Hermitian
    if isfield(Array, 'alpha_types'), Array.alpha_types = real(Array.alpha_types);
    else,                             Array.alpha = real(Array.alpha);
    end
end
% One polarizability per particle type (dda_alpha reads it); the old solvers get
% the full array back with dda_expand_alpha.
Array = dda_compact_alpha(Array);
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

function Array = mie_per_radius(Array, Input)
% Mie polarizability for each distinct radius in Array.radius (array editor)
[ru, ~, typ] = unique(Array.radius(:));
types = [];
for m = 1:numel(ru)
    In1 = Input;  In1.particle.r = ru(m);
    A1 = mie_one_type(Array, In1);
    if m == 1, Array = A1; types = A1.alpha_types; else, types = cat(1, types, A1.alpha_types); end
end
Array.alpha_types = types;
Array.alpha_type  = typ(:);
end

function Array = mie_one_type(Array, Input)
% Mie polarizability of one particle, shared by all (identical particles)
In1 = Input;   In1.N = 1;
A1  = Array;   A1.r_i = Array.r_i(1,:);
A1  = Mie_polarizability_DDA(A1, In1);
f = fieldnames(A1);
for n = 1:numel(f)
    if ~any(strcmp(f{n}, {'r_i', 'alpha'})), Array.(f{n}) = A1.(f{n}); end
end
Array.alpha_types = A1.alpha;                 % 1 x 3 x 3 x Nk (or 1 x Nk for cyl_TM)
Array.alpha_type  = ones(size(Array.r_i,1), 1);
if isfield(Array, 'alpha'), Array = rmfield(Array, 'alpha'); end
end
