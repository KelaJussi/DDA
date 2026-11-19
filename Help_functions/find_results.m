function varargout = find_results(varargin)
% FIND_RESULTS  Search saved results by their Input parameters.
%
%   T = find_results()                       % all complex-frequency results, newest first
%   T = find_results('Kind','real', ...)     % real-frequency results
%   T = find_results('Kind','all', ...)      % both
%   T = find_results('Name', value, ...)     % only matching results
%   find_results(...)                        % print a numbered list instead
%   plot_results(T, 1)                       % plot row 1 with the right plot function
%
% Every table has the columns kind ("complex" or "real"), folder, file, date
% and bytes, then one column per Input field. Real-frequency results also have
% info_* columns (info_created, info_solver, info_run_solve_time_s, ...).
%
% Names are Input fields; nested fields with '.' or '_' ('particle.r' or
% 'particle_r'). Long arrays such as range.real_parts are indexed as
% range_real_parts_N, _min and _max (and _imin, _imax if complex).
% Columns 'file' and 'date' can be filtered too.
%
% Value types:
%   number / vector  equal within relative tolerance Tol (default 1e-6),
%                    vectors must have the same length:  'p', [300e-9 300e-9]
%   text             case-insensitive, * and ? wildcards: 'save_message', '*test*'
%   cell array       any of the values matches:          'dimensionality', {'2D','3D'}
%   function handle  custom test on the stored value:    'particle_r', @(r) r >= 40e-9
%                                                        'date', @(d) d > datetime(2026,10,1)
% Options (same name/value form):
%   'Kind'     'complex' (default), 'real' or 'all'
%   'Folder'   results folder (default dda_results_folder(Kind); not with 'all')
%   'Rebuild'  true = re-read every file instead of using the cache
%   'Tol'      relative tolerance for numbers
%   'KeepEmpty' true = keep columns that are empty for all matches
%
% Files saved before an option existed get the values listed in
% "defaults" in build_results_index (e.g. periodic_boundary = [0 0]).
%
% Example:
%   T = find_results('dimensionality','3D', 'arrayN',[20 20], 'periodic_boundary',[1 0]);
%   T(:, {'file','date','p','particle_r'})
%   plot_results(T, 1)
%
% See README_results_index.md.

%% Parse arguments
opts = struct('Kind', 'complex', 'Folder', '', 'Rebuild', false, 'Tol', 1e-6, 'KeepEmpty', false);
optnames = fieldnames(opts);
if mod(numel(varargin), 2) ~= 0
    error('find_results:args', 'Arguments must be name/value pairs.');
end
filters = cell(0, 2);
for k = 1:2:numel(varargin)
    name = char(varargin{k});
    o = find(strcmpi(name, optnames), 1);
    if ~isempty(o)
        opts.(optnames{o}) = varargin{k+1};
    else
        filters(end+1, :) = {name, varargin{k+1}}; %#ok<AGROW>
    end
end
kind = lower(string(opts.Kind));
if ~ismember(kind, ["complex","real","all"])
    error('find_results:kind', 'Kind must be ''complex'', ''real'' or ''all''.');
end

%% Index and filter
if kind == "all"
    if ~isempty(opts.Folder), error('find_results:folder', 'Folder cannot be combined with Kind ''all''.'); end
    cols = stack_columns(index_of("complex", '', opts.Rebuild), index_of("real", '', opts.Rebuild));
    opts.Folder = 'all result folders';
else
    cols = index_of(kind, opts.Folder, opts.Rebuild);
    if isempty(opts.Folder), opts.Folder = dda_results_folder(kind); end
end
ntot = numel(cols.file);
datenums = cols.date;
try
    cols.date = datetime(cols.date, 'ConvertFrom', 'datenum');
catch
    % no datetime available: keep datenum
end

names = fieldnames(cols);
mask = true(ntot, 1);
for f = 1:size(filters, 1)
    col = resolve_column(filters{f,1}, names);
    mask = mask & match_column(cols.(col), filters{f,2}, opts.Tol);
end

% Keep matches, newest first
idx = find(mask);
[~, ord] = sort(datenums(idx), 'descend');
idx = idx(ord);
datenums = datenums(idx);
for k = 1:numel(names)
    c = cols.(names{k})(idx);
    if ~opts.KeepEmpty && ~any(strcmp(names{k}, {'kind','folder','file','date','bytes'})) && ~isempty(idx) && all_missing(c)
        cols = rmfield(cols, names{k});
    else
        cols.(names{k}) = c;
    end
end

%% Output
if nargout > 0
    T = struct2table(cols);
    T.Properties.UserData = struct('folder', opts.Folder);
    varargout{1} = T;
else
    fprintf('%d of %d results match (folder: %s)\n', numel(idx), ntot, opts.Folder);
    for i = 1:numel(idx)
        fprintf('%4d  %-7s %s  %s\n', i, cols.kind{i}, datestr(datenums(i), 'yyyy-mm-dd HH:MM'), cols.file{i});
    end
end
end

%% ------------------------------------------------------------------------
function cols = index_of(kind, folder, rebuild)
% Index of one result folder, with kind and folder columns first
if isempty(folder), folder = dda_results_folder(kind); end
c = build_results_index(folder, rebuild);
n = numel(c.file);
cols = struct('kind', {repmat({char(kind)}, n, 1)}, 'folder', {repmat({char(folder)}, n, 1)});
f = fieldnames(c);
for k = 1:numel(f), cols.(f{k}) = c.(f{k}); end
end

function c = stack_columns(a, b)
% Put two indexes below each other; a column missing in one gets NaN / '' / []
na = numel(a.file); nb = numel(b.file);
names = unique([fieldnames(a); fieldnames(b)], 'stable');
c = struct();
for k = 1:numel(names)
    nm = names{k};
    ca = get_col(a, nm, na, b); cb = get_col(b, nm, nb, a);
    if isnumeric(ca) && isnumeric(cb)
        c.(nm) = [ca; cb];
    else
        if ~iscell(ca), ca = num2cell(ca); end
        if ~iscell(cb), cb = num2cell(cb); end
        c.(nm) = [ca; cb];
    end
end
end

function col = get_col(s, nm, n, other)
if isfield(s, nm)
    col = s.(nm);
elseif isnumeric(other.(nm))
    col = nan(n, 1);
elseif iscell(other.(nm)) && ~isempty(other.(nm)) && all(cellfun(@ischar, other.(nm)))
    col = repmat({''}, n, 1);
else
    col = cell(n, 1);
end
end

function col = resolve_column(name, names)
key = strrep(char(name), '.', '_');
if strncmp(key, 'Input_', 6) && ~any(strcmp(key, names))
    key = key(7:end);                      % allow 'Input.p'
end
i = find(strcmp(key, names), 1);
if isempty(i), i = find(strcmpi(key, names), 1); end
if ~isempty(i)
    col = names{i};
    return
end
similar = names(~cellfun(@isempty, strfind(lower(names), lower(key))));
if isempty(similar)
    hint = 'See T = find_results(); T.Properties.VariableNames for all columns.';
else
    hint = ['Did you mean: ' strjoin(similar(:).', ', ') '?'];
end
error('find_results:unknownField', 'No column "%s" in the results index. %s', key, hint);
end

function m = match_column(col, q, tol)
n = size(col, 1);
m = false(n, 1);
for i = 1:n
    if iscell(col), v = col{i}; else, v = col(i); end
    m(i) = match_value(v, q, tol);
end
end

function m = match_value(v, q, tol)
m = false;
if iscell(q)                                   % any of the alternatives
    for j = 1:numel(q)
        if match_value(v, q{j}, tol), m = true; return, end
    end
elseif isa(q, 'function_handle')
    if isempty(v) || (isnumeric(v) && isscalar(v) && isnan(v)), return, end
    try
        r = q(v);
        m = ~isempty(r) && all(logical(r(:)));
    catch
        m = false;
    end
elseif ischar(q) || isstring(q)
    if ~ischar(v) || isempty(v), return, end
    q = char(q);
    if any(q == '*' | q == '?')
        m = ~isempty(regexpi(v, ['^' regexptranslate('wildcard', q) '$'], 'once'));
    else
        m = strcmpi(v, q);
    end
elseif isnumeric(q) || islogical(q)
    if ~(isnumeric(v) || islogical(v)) || numel(v) ~= numel(q) || isempty(v), return, end
    v = double(v(:)); q = double(q(:));
    d = abs(v - q);
    m = all(d == 0 | d <= tol*max(abs(v), abs(q)));
end
end

function tf = all_missing(c)
if iscell(c)
    tf = all(cellfun(@isempty, c));
elseif isnumeric(c)
    tf = all(isnan(c));
else
    tf = false;
end
end
