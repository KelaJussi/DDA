function varargout = dda_favorites(action, varargin)
% DDA_FAVORITES  Saved input profiles ("favorites") for DDA runs.
%
%   folder = dda_favorites('folder')                 where favorites are kept
%   list   = dda_favorites('list')                   struct array: name, note, saved, n, layout, file
%   dda_favorites('save', name, note, paths, values, array_edit)
%   fav    = dda_favorites('load', name)
%   dda_favorites('delete', name)
%   Input  = dda_favorites('apply', Input, name)     put a favorite's values into Input
%
% A favorite stores the value of every option of the parameter GUI (paths such
% as "particle.r" with their values, in SI units) and, if present, the edited
% array layout (Input.array_edit). Loading it reproduces those settings exactly,
% also if the defaults in Make_predefined_Inputs change later.
%
% Favorites are .mat files in DDA_favorites next to DDA_Master.m, one per
% favorite. Use them from the Favorites tab of the GUI, or in a script with
% Input.favorite = "name" in Make_predefined_Inputs.

switch lower(string(action))
    case "folder"
        varargout{1} = fav_folder();

    case "list"
        d = dir(fullfile(fav_folder(), '*.mat'));
        list = struct('name', {}, 'note', {}, 'saved', {}, 'n', {}, 'layout', {}, 'file', {});
        for k = 1:numel(d)
            try
                S = load(fullfile(d(k).folder, d(k).name), 'fav');
                f = S.fav;
                list(end+1) = struct('name', string(f.name), 'note', string(f.note), ...
                    'saved', f.saved, 'n', numel(f.paths), 'layout', ~isempty(f.array_edit), ...
                    'file', string(fullfile(d(k).folder, d(k).name))); %#ok<AGROW>
            catch
                warning('dda_favorites:unreadable', 'Skipping unreadable favorite %s', d(k).name);
            end
        end
        if ~isempty(list)
            [~, o] = sort(lower([list.name]));  list = list(o);
        end
        varargout{1} = list;

    case "save"
        [name, note, paths, values, array_edit] = varargin{1:5};
        name = strtrim(string(name));
        if strlength(name) == 0, error('dda_favorites:name', 'Give the favorite a name.'); end
        f = fav_file(name);
        if isfile(f)
            old = load(f, 'fav');
            if lower(string(old.fav.name)) ~= lower(name)
                error('dda_favorites:clash', ['"%s" would use the same file as the favorite "%s" ' ...
                    '(names differing only in special characters). Choose another name.'], name, old.fav.name);
            end
        end
        fav = struct('name', name, 'note', string(note), 'saved', datetime('now'), ...
                     'paths', string(paths(:)), 'values', {values(:)}, ...
                     'array_edit', array_edit, 'version', 1);
        save(fav_file(name), 'fav');
        varargout{1} = fav_file(name);

    case "load"
        f = fav_file(varargin{1});
        if ~isfile(f), error('dda_favorites:missing', 'There is no favorite "%s" (looked for %s).', varargin{1}, f); end
        S = load(f, 'fav');
        varargout{1} = S.fav;

    case "delete"
        f = fav_file(varargin{1});
        if isfile(f), delete(f); end

    case "apply"
        Input = varargin{1};
        fav = dda_favorites('load', varargin{2});
        for k = 1:numel(fav.paths)
            parts = cellstr(split(fav.paths(k), "."));
            Input = setfield(Input, parts{:}, fav.values{k});
        end
        if isempty(fav.array_edit)
            if isfield(Input, 'array_edit'), Input = rmfield(Input, 'array_edit'); end
        else
            Input.array_edit = fav.array_edit;
        end
        fprintf('Favorite "%s" applied (%d setting(s)%s).\n', fav.name, numel(fav.paths), ...
            ternary(isempty(fav.array_edit), "", " and an edited array layout"));
        varargout{1} = Input;

    otherwise
        error('dda_favorites:action', 'Unknown action "%s".', action);
end
end

%% ------------------------------------------------------------------------
function f = fav_folder()
root = fileparts(which('DDA_Master'));
if isempty(root), root = pwd; end
f = fullfile(root, 'DDA_favorites');
if ~isfolder(f), mkdir(f); end
end

function f = fav_file(name)
% file name from the favorite's name: letters, digits, '-' and '_' only
s = regexprep(char(strtrim(string(name))), '[^A-Za-z0-9_\-]', '_');
f = fullfile(fav_folder(), [s '.mat']);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
