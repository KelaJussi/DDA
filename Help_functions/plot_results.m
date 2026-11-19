function plot_results(T, rows, plotting)
% PLOT_RESULTS  Plot results found with find_results.
%
%   T = find_results('Kind', 'all', 'arrayN', [20 20 1]);
%   plot_results(T)              % row 1
%   plot_results(T, [1 3])       % rows 1 and 3
%   plot_results(T, 2, opts)     % with plot-option overrides:
%                                %   real:    a plotting struct, e.g. struct('Q', struct('enable', false))
%                                %   complex: complex options, e.g. struct('quiver', false)
%
% Chooses real_freqs_plots_from_file or complex_freqs_plots_from_file from the
% kind column of T, and opens the file in the folder it was found in.

if nargin < 2 || isempty(rows), rows = 1; end
if nargin < 3, plotting = []; end
names = T.Properties.VariableNames;
for r = rows(:).'
    if r < 1 || r > height(T)
        warning('plot_results:row', 'There is no row %d (the table has %d).', r, height(T));
        continue
    end
    kind = "complex";
    if ismember('kind', names), kind = string(T.kind{r}); end
    if ismember('folder', names)
        file = fullfile(T.folder{r}, T.file{r});
    else
        file = fullfile(dda_results_folder(kind), T.file{r});
    end
    fprintf('Plotting %s result: %s\n', kind, T.file{r});
    if kind == "real"
        real_freqs_plots_from_file(file, plotting);
    else
        complex_freqs_plots_from_file(file, plotting);
    end
end
end
