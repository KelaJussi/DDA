function folder = dda_results_folder(kind)
% DDA_RESULTS_FOLDER  Folder where results are saved.
%
%   dda_results_folder()            complex-frequency results  (DDA/Complex_data)
%   dda_results_folder("complex")   same
%   dda_results_folder("real")      real-frequency results     (DDA/Real_freq_results)
%
% The DDA main folder is the parent of Help_functions (where this file lives),
% so this works wherever the code is. Both folders are created on first use
% (and are git-ignored).
if nargin < 1, kind = "complex"; end
root = fileparts(fileparts(mfilename('fullpath')));

switch lower(string(kind))
    case "complex"
        folder = fullfile(root, 'Complex_data');
        if ~isfolder(folder), mkdir(folder); end
    case "real"
        folder = fullfile(root, 'Real_freq_results');
        if ~isfolder(folder), mkdir(folder); end
    otherwise
        error('dda_results_folder:kind', 'kind must be "complex" or "real".');
end
end
