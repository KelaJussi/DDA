function filename = save_complex_freqs_modes(C,Input,Array, filename)
    if nargin < 4 || isempty(filename)   % was nargin < 3, so the name was never constructed
        filename = construct_filename(Input);
    end
    folder = dda_results_folder();
    if isfield(C,'P_eigs')
        C = rmfield(C,'P_eigs');
    end
    if isfield(C,'vectors_all')
        C = rmfield(C,'vectors_all');
    end
    if isfield(Array,'alpha')
        Array = rmfield(Array,'alpha');
    end
    if isfield(Array,'gamma')
        Array = rmfield(Array,'gamma');
    end
    fullpath = fullfile(folder, filename);
    save(fullpath,"C","Input","Array")
    disp(append("saved results to ", fullpath))
    % The results index (find_results) picks up the new file automatically.
end
