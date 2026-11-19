function [P,C,Field,Array] = Calculate_DDA_eigenmodes(Input,Array)
% CALCULATE_DDA_EIGENMODES  One eigenmode solver for every dimensionality.
%
% Replaces:
%   Calculate_1D_DDA_complex_freqs_eigenmodes   Calculate_1D_DDA_generic_eigenmodes
%   Calculate_2D_DDA_complex_freqs_eigenmodes   Calculate_2D_DDA_generic_eigenmodes
%   Calculate_3D_DDA_complex_freqs_eigenmodes   Calculate_3D_DDA_generic_eigenmodes
%   Calculate_cyl_TM_DDA_complex_freqs_eigenmodes (and the missing cyl_TM generic one)
%
% Real-frequency runs are simply Input.range.imag_parts = 0.
%
% Options (all optional, defaults in brackets; set them in DDA_Master):
%   Input.parallel        "none" | "threads"
%                         [legacy: "threads" if Input.useparallel, else "none"]
%                         Cluster and background runs are submitted by
%                         submit_dda_job, which runs this solver inside the job.
%   Input.threads.*       num_workers [] (= all cores), keep_pool [true]
%   Input.track.*         see Track_eigenmodes.m; enable ["auto" = real-frequency runs only]
%   Input.degeneracy_tol  [1e29] absolute tolerance used for C.degeneracy
%   Input.P_method        ["direct"] how the driven response P is computed:
%                         "direct"     P = A\Einc, exact and robust also near
%                                      degenerate modes (ill-conditioned eigenvectors)
%                         "eigenmodes" P = V*(c./lambda), the sum over eigenmodes
%                                      (same result when V is well conditioned;
%                                      shows how the eigenmodes build the response)
%
% Outputs
%   P      [Ndof x 1 x Nk] single, driven polarization (Input.P_method)
%   C      eigenvalues, smallest eigenvalue/vector, projections of Einc,
%          extinction, and (if tracking) mode-ordered fields for JTK_plots
%   Field  incident field at the last frequency (as before), plus
%          Field.Einc_all [Ndof x Nk], the incident field at every frequency
%
% DRAFT - written without MATLAB at hand. Validate against the old solvers
% (Input.use_legacy_eig_solvers = true in DDA_Master) before deleting them.

ncomp = dda_ncomp(Input);
Ndof  = ncomp*Input.N;
Nk    = Input.range.k_N;
track = track_enabled(Input);

sz   = check_sizes(Input);
lean = sz.lean;      % complex frequencies: keep only what complex_freqs_plots_from_file uses

switch parallel_mode(Input)
    case "threads"
        topt = struct_defaults(getopt(Input,'threads',struct()), ...
                               struct('num_workers',[],'keep_pool',true));
        pool = start_thread_pool(topt.num_workers);
        [C,P] = run_core(Input,Array,Ndof,track,lean,pool.NumWorkers);
        if ~topt.keep_pool, delete(pool); end
    otherwise
        [C,P] = run_core(Input,Array,Ndof,track,lean,0);   % parfor with M = 0: serial, never starts a pool
end

Field = Make_3D_Field_wl_looped(Input,Array,Nk);
Field.Einc_all = C.Einc;            % Ndof x Nk, saved with real-frequency results ([] for lean runs)
C = rmfield(C,'Einc');
end

%% ========================================================================
%  Core
%% ========================================================================
function [C,P] = run_core(Input,Array,Ndof,track,lean,maxWorkers)
S = solve_all_k(Input,Array,Ndof,track,lean,maxWorkers);

C = struct();
C.solver_version  = 4;          % 2: matrix-inverse self term (see Make_int_matrix_A_SI)
                                % 3: exact Mie polarizability of the cylinders (cyl_TM)
                                % 4: mirrors no longer double-count particles on the mirror plane;
                                %    Make_int_matrix_A_SI in double precision (was single)
if lean
    % Complex frequencies: the smallest eigenvalue, its eigenvector and the
    % degeneracy are all complex_freqs_plots_from_file uses. The other
    % eigenvalues, the coefficients, P and the incident field are not computed
    % (Input.save.complex_all = true keeps them).
    C.min_eig_values  = S.min_vals;
    C.min_eig_vectors = S.min_vecs;
    C.degeneracy      = S.degeneracy;
    C.lean            = true;
    C.Einc            = [];
    P = [];
    return
end
C.values          = S.values;
C.min_eig_values  = S.min_vals;
C.min_eig_vectors = S.min_vecs;
C.degeneracy      = S.degeneracy;
C.coefficients    = S.coefs;
[~, C.fundindex]  = max(abs(S.coefs),[],1);

C = Compute_mode_extinction(C, Input, Array);
C.ext  = dda_extinction(S.P, S.Einc, Input, Array);   % exact: P = A\Einc
C.Einc = single(S.Einc);

if track
    % Done here (not on the client) so a cluster job only sends back results,
    % not every eigenvector of every frequency.
    C = Track_eigenmodes(C, S.vectors, Input, Array);
end

P = reshape(single(S.P), Ndof, 1, []);
end

function S = solve_all_k(Input,Array,Ndof,keep_vectors,lean,maxWorkers)
Nk  = Input.range.k_N;
tol = getopt(Input,'degeneracy_tol',1e29);
P_direct = ~strcmpi(string(getopt(Input,'P_method',"direct")), "eigenmodes");

% M must be a nonnegative integer. Inf means "use the pool that is already
% open"; with no open pool fall back to 0 (serial), so this function never
% triggers parfor's automatic pool creation.
if isinf(maxWorkers)
    pool = gcp('nocreate');
    if isempty(pool), maxWorkers = 0; else, maxWorkers = pool.NumWorkers; end
end

nfull      = Ndof*~lean;           % lean runs keep no per-eigenvector arrays
values     = complex(NaN(nfull,Nk));
min_vecs   = complex(NaN(Ndof,Nk,'single'));
coefs      = complex(NaN(nfull,Nk));
Pk         = complex(NaN(nfull,Nk));
Ek         = complex(NaN(nfull,Nk));
min_vals   = complex(NaN(1,Nk));
degeneracy = NaN(1,Nk);
vectors    = cell(1,Nk);
t0 = tic;

parfor (kk = 1:Nk, maxWorkers)
    A      = Make_int_matrix(Input,Array,kk);
    [V,D]  = eig(A);
    d      = diag(D);
    [mv,i] = min(d);                       % complex min compares |d|

    min_vals(kk)   = mv;
    min_vecs(:,kk) = single(V(:,i));
    degeneracy(kk) = nnz(abs(d - mv) <= tol);
    if ~lean
        Fk = Make_3D_Field_wl_looped(Input,Array,kk);
        ck = V \ double(Fk.Einc(:));
        values(:,kk) = d;
        coefs(:,kk)  = ck;
        if P_direct
            Pk(:,kk) = A \ double(Fk.Einc(:));   % driven response, exact
        else
            Pk(:,kk) = V*(ck./d);             % the same as a sum over eigenmodes
        end
        Ek(:,kk)     = double(Fk.Einc(:));
    end
    if keep_vectors, vectors{kk} = V; end

    Display_calc_times_and_convergence_specs('fwm100',kk,Array.k,' A inverted in ',0,NaN,t0);
end

S = struct('values',values,'min_vals',min_vals,'min_vecs',min_vecs, ...
           'degeneracy',degeneracy,'coefs',coefs,'P',Pk,'Einc',Ek);
S.vectors = vectors;
end

%% ========================================================================
%  Local helpers
%% ========================================================================
function pool = start_thread_pool(num_workers)
pool = gcp('nocreate');
if isempty(pool)
    if isempty(num_workers), pool = parpool("Threads");
    else,                    pool = parpool("Threads", num_workers);
    end
end
end

function C = Compute_mode_extinction(C, Input, Array)
% Extinction of each eigenmode on its own: what dda_extinction gives when the
% incident field is that (unit-norm) eigenvector, i.e. 4*pi*k*Im(1/lambda)/cell.
% Same normalisation as dda_extinction, so C.extsall and C.ext are comparable.
eps0 = 8.8541878188e-12;
if any(Input.dimensionality == ["cyl_TM","cyl_TE"]), cell_area = Input.p(1);
else,                                              cell_area = Input.p(1)*Input.p(2);
end
exts = double(4*pi.*Array.k(:).'.*imag(1./C.values)./cell_area);
if dda_alpha_is_SI(Input, Array)
    exts = exts/4/pi/eps0./double(Array.n(:).').^2;
end
C.exts_unordered = exts;
end

function A = Make_int_matrix(Input,Array,kk)
if Input.dimensionality == "cyl_TM"
    A = Make_cyl_TM_matrix_A_SI(Input,Array,kk);   % + layered background (Input.layers)
elseif Input.dimensionality == "cyl_TE"
    A = Make_cyl_TE_matrix_A_SI(Input,Array,kk);   % + layered background (Input.layers)
else
    A = Make_int_matrix_A_SI(Input,Array,kk);   % handles 1Dx/1Dy/1Dz/2D/3D
end
end

function n = dda_ncomp(Input)
switch string(Input.dimensionality)
    case "3D", n = 3;
    case {"2D","cyl_TE"}, n = 2;
    otherwise, n = 1;                            % 1Dx, 1Dy, 1Dz, cyl_TM
end
end

function m = parallel_mode(Input)
if isfield(Input,'parallel'), m = string(Input.parallel);
elseif getopt(Input,'useparallel',false), m = "threads";
else, m = "none";
end
if ~ismember(m, ["none","threads"])
    error('Calculate_DDA_eigenmodes:parallel', ...
        ['Input.parallel = "%s" is handled by submit_dda_job (via DDA_calc_looped); ' ...
         'this solver itself runs "none" or "threads".'], m);
end
end

function tf = track_enabled(Input)
e = getopt(getopt(Input,'track',struct()),'enable',"auto");
if isstring(e) || ischar(e)
    tf = ~any(Input.range.imag_parts);           % "auto": real-frequency runs only
else
    tf = logical(e);
end
end

function sz = check_sizes(Input)
% Same estimate as dda_size_check, which runs before jobs are submitted
sz = dda_result_size(Input);
if sz.track
    gb  = sz.vectors/1e9;                        % all eigenvectors are kept for tracking
    lim = getopt(getopt(Input,'track',struct()),'max_memory_GB',8);
    if gb > lim
        error(['Mode tracking needs %.1f GB for the eigenvectors (limit %.1f GB). ' ...
               'Reduce the frequency grid, set Input.track.enable = false, or raise ' ...
               'Input.track.max_memory_GB.'], gb, lim);
    end
end
if any(Input.range.imag_parts) && sz.results > 1.8e9   % saved by save_complex_freqs_modes (v7 format)
    error(['Complex-frequency results would be %.2f GB, over the 2 GB per-variable limit of ' ...
           '.mat files. Split the frequency grid into several runs.'], sz.results/1e9);
end
end

function s = struct_defaults(s, d)
f = fieldnames(d);
for n = 1:numel(f)
    if ~isfield(s,f{n}), s.(f{n}) = d.(f{n}); end
end
end

function v = getopt(s,name,default)
if isstruct(s) && isfield(s,name), v = s.(name); else, v = default; end
end
