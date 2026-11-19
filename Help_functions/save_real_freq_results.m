function fullpath = save_real_freq_results(fullpath, C, P, Field, Input, Array)
% SAVE_REAL_FREQ_RESULTS  Store a real-frequency run so real_freqs_plots_from_file
% can redo every plot without recalculating.
%
% Saved variables
%   Input  the full Input struct (incl. Input.plotting at run time)
%   C      solver results: ext, and for the eigen solver values,
%          coefficients, extsall, coefs_ordered, P_eigs, track ...
%   P      Ndof x Nk dipole moments (driven response A\Einc)
%   Einc   Ndof x Nk incident field at every frequency
%   Array  positions, k, n, wl, alpha (no function handles / GPU arrays)
%   info   date, solver, MATLAB version, file format version
%
% Options (Input.save.*, defaults in brackets)
%   P_eigs       [true] keep per-mode dipole moments (Ndof x Nk per mode)
%   max_file_GB  [4]    drop P_eigs, with a warning, if the file would be larger
%
% Saved with -v7.3, so single variables are not limited to 2 GB.

opt = struct('P_eigs', true, 'max_file_GB', 4);
if isfield(Input,'save'), opt = merge(opt, Input.save); end

Nk = Input.range.k_N;
P  = reshape(P, [], Nk);

%% Incident field at every frequency
if isfield(Field,'Einc_all')
    Einc = Field.Einc_all;
else                                  % e.g. Calculate_3D_DDA_generic returns one frequency
    Einc = complex(zeros(size(P),'single'));
    for kk = 1:Nk
        F = Make_3D_Field_wl_looped(Input, Array, kk);
        Einc(:,kk) = F.Einc(:);
    end
end

%% One extinction formula for every solver
% (not for the legacy eigen solvers: their P is the fundamental mode, not A\Einc)
legacy = isfield(Input,'use_legacy_eig_solvers') && Input.use_legacy_eig_solvers ...
         && string(Input.inversion_tool) == "eig";
if ~legacy
    if string(Input.inversion_tool) ~= "eig" && isfield(C,'ext')
        C.ext_solver = C.ext;         % keep the direct solver's own value
    end
    C.ext = dda_extinction(P, Einc, Input, Array);
end

%% Drop what cannot or need not be saved
Array = strip_unsaveable(Array);
C     = strip_unsaveable(C);
if isfield(C,'P_out'),     C = rmfield(C,'P_out');     end   % = P, saved separately
if isfield(C,'Field_out'), C = rmfield(C,'Field_out'); end

if isfield(C,'P_eigs')
    w = whos('C');
    if ~opt.P_eigs || w.bytes/1e9 > opt.max_file_GB
        if opt.P_eigs
            warning('save_real_freq_results:size', ...
                ['Results are %.1f GB; C.P_eigs is not saved (Input.save.max_file_GB = %g). ' ...
                 'Mode dipole plots will be skipped.'], w.bytes/1e9, opt.max_file_GB);
        end
        C = rmfield(C,'P_eigs');
    end
end

info = struct('created', datetime('now'), 'solver', string(Input.inversion_tool), ...
              'matlab', version, 'format', 1);
if isfield(C, 'run_info')
    info.run = C.run_info;      % solve time, workers, threads, node ... (indexed by find_results)
end

P    = single(P);     %#ok<NASGU>
Einc = single(Einc);  %#ok<NASGU>
if strlength(string(fullpath)) > 250
    warning('save_real_freq_results:longPath', ...
        'Path is %d characters; Windows may refuse paths over 260.', strlength(string(fullpath)));
end
save(fullpath, 'Input', 'C', 'P', 'Einc', 'Array', 'info', '-v7.3');
fprintf('Saved real-frequency results to %s\n', fullpath);
end

%% ------------------------------------------------------------------------
function s = strip_unsaveable(s)
f = fieldnames(s);
for n = 1:numel(f)
    v = s.(f{n});
    if isa(v,'function_handle')
        s = rmfield(s, f{n});
    elseif isa(v,'gpuArray')
        s.(f{n}) = gather(v);
    end
end
end

function a = merge(a, b)
f = fieldnames(b);
for n = 1:numel(f), a.(f{n}) = b.(f{n}); end
end
