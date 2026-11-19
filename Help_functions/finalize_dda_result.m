function finalize_dda_result(R, info, plot_now)
% FINALIZE_DDA_RESULT  Save a finished DDA run and optionally plot it.
% The same step for blocking runs (DDA_calc_looped) and for collected jobs
% (collect_dda_jobs).
%
%   R     output of dda_solve_job: R.P, R.C, R.Field, R.Array
%   info  struct with fullpath, filename, kind ("real"/"complex"), Input
%         (a registry entry has all of these)

if isfield(R, 'info')
    R.C.run_info = R.info;      % saved with the results: workers, threads, times
end

switch string(info.kind)
    case "real"
        save_real_freq_results(char(info.fullpath), R.C, R.P, R.Field, info.Input, R.Array);
    case "complex"
        save_complex_freqs_modes(R.C, info.Input, R.Array, info.filename);
    otherwise
        error('finalize_dda_result:kind', 'Unknown result kind "%s".', info.kind);
end

if plot_now
    if string(info.kind) == "real"
        real_freqs_plots_from_file(char(info.fullpath));
    else
        complex_freqs_plots_from_file(info.filename);
    end
end
end
