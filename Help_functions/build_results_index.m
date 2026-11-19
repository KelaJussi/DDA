function cols = build_results_index(folder, rebuild)
% BUILD_RESULTS_INDEX  Index of the Input parameters of all saved results.
%
%   cols = build_results_index()            % default folder, incremental
%   cols = build_results_index(folder, true) % re-read every file
%
% Reads only the Input variable (and, in real-frequency files, the small info
% variable with creation date and run information) from each .mat file in the
% folder and flattens it to one row per file (nested fields become e.g.
% particle_r; info fields become info_created, info_run_solve_time_s, ...).
% Works for both result folders: dda_results_folder("complex") (default) and
% dda_results_folder("real"); find_results chooses the folder.
% The rows are cached in <folder>\results_index.mat; on later calls only
% new or changed files are read and deleted files are dropped, so the
% index always matches what is on disk.
%
% Returns a struct of equal-length columns (find_results turns it into a
% table). Normally you call find_results instead of this function.

if nargin < 1 || isempty(folder), folder = dda_results_folder(); end
if nargin < 2, rebuild = false; end
folder = char(folder);

INDEX_VERSION = 2;   % bump if the flattening rules below change (2: info, datetimes)
cachefile = fullfile(folder, 'results_index.mat');

% Values assumed for files saved before an option existed.
% Add an entry here whenever you add a new Input option.
defaults = struct( ...
    'periodic_boundary', [0 0], ...
    'proj_real_freq',    0);

%% Load cache
entries = empty_entries();
if ~rebuild && isfile(cachefile)
    try
        S = load(cachefile, 'entries', 'index_version');
        if isfield(S,'index_version') && S.index_version == INDEX_VERSION
            entries = S.entries;
        end
    catch
        warning('build_results_index:cache', 'Could not read %s, rebuilding the index.', cachefile);
    end
end

%% Compare folder contents with cache
d = dir(fullfile(folder, '*.mat'));
d = d(~[d.isdir] & ~strcmp({d.name}, 'results_index.mat'));
n = numel(d);

cached_names = {entries.name};
src = zeros(n,1);   % index into cache, 0 = (re)read
for i = 1:n
    j = find(strcmp(cached_names, d(i).name), 1);
    if ~isempty(j) && entries(j).datenum == d(i).datenum && entries(j).bytes == d(i).bytes
        src(i) = j;
    end
end

nread = nnz(src == 0);
if nread > 0
    fprintf('Indexing %d new or changed result file(s) in %s ...\n', nread, folder);
end

newlist = cell(n,1);
for i = 1:n
    if src(i) > 0
        newlist{i} = entries(src(i));
    else
        newlist{i} = read_entry(folder, d(i));
    end
end
if n > 0
    new_entries = [newlist{:}];
else
    new_entries = empty_entries();
end

%% Save cache if anything changed
if nread > 0 || numel(entries) ~= n
    entries = new_entries; index_version = INDEX_VERSION; %#ok<NASGU>
    try
        save(cachefile, 'entries', 'index_version');
    catch err
        warning('build_results_index:save', 'Could not save index cache: %s', err.message);
    end
end
entries = new_entries;

% Report files without a usable Input (only when they are first read)
justread = entries(src == 0);
bad = justread(~[justread.ok]);
for i = 1:numel(bad)
    warning('build_results_index:skipped', 'Not indexed: %s (%s)', bad(i).name, bad(i).msg);
end

cols = merge_entries(entries([entries.ok]), defaults);
end

%% ------------------------------------------------------------------------
function e = empty_entries()
e = struct('name',{}, 'datenum',{}, 'bytes',{}, 'ok',{}, 'vals',{}, 'msg',{});
end

function e = read_entry(folder, f)
e = struct('name', f.name, 'datenum', f.datenum, 'bytes', f.bytes, ...
           'ok', false, 'vals', struct(), 'msg', '');
w = warning('off', 'all');            % silence "variable not found" etc.
restore = onCleanup(@() warning(w));
try
    vars = who('-file', fullfile(folder, f.name));
    want = intersect({'Input','info'}, vars);
    S = load(fullfile(folder, f.name), want{:});   % only Input (and info), not the results
    if ~isfield(S, 'Input') || ~isstruct(S.Input) || ~isscalar(S.Input)
        e.msg = 'no Input variable';
        return
    end
    e.vals = flatten_into(struct(), S.Input, '');
    if isfield(S, 'info') && isstruct(S.info) && isscalar(S.info)
        e.vals = flatten_into(e.vals, S.info, 'info_');
    end
    e.ok = true;
catch err
    e.msg = err.message;
end
end

function vals = flatten_into(vals, s, prefix)
% Flatten a (nested) struct into one level of scalar / short-vector values.
MAXVEC = 10;   % longer numeric arrays are summarised as _N, _min, _max
fn = fieldnames(s);
for i = 1:numel(fn)
    v = s.(fn{i});
    name = [prefix fn{i}];
    if isstruct(v) && isscalar(v)
        vals = flatten_into(vals, v, [name '_']);
    elseif ischar(v) || isstring(v)
        vals = setval(vals, name, strjoin(cellstr(v), ', '));
    elseif isdatetime(v) && isscalar(v)
        if isnat(v), vals = setval(vals, name, ''); else, vals = setval(vals, name, char(v, 'yyyy-MM-dd HH:mm:ss')); end
    elseif (isnumeric(v) || islogical(v)) && ~isempty(v)
        v = double(v);
        if numel(v) <= MAXVEC
            vals = setval(vals, name, reshape(v, 1, []));
        else
            r = real(v(:));
            vals = setval(vals, [name '_N'],   numel(v));
            vals = setval(vals, [name '_min'], min(r));
            vals = setval(vals, [name '_max'], max(r));
            if ~isreal(v)
                vals = setval(vals, [name '_imin'], min(imag(v(:))));
                vals = setval(vals, [name '_imax'], max(imag(v(:))));
            end
        end
    end
    % cells, function handles and objects are not indexed
end
end

function vals = setval(vals, name, v)
if any(strcmp(name, {'file','date','bytes'}))
    name = ['Input_' name];           % avoid clashing with the file columns
end
name = name(1:min(end, namelengthmax));
vals.(name) = v;
end

function cols = merge_entries(entries, defaults)
% Combine rows with possibly different fields into equal-length columns.
n = numel(entries);
cols = struct();
cols.file  = reshape({entries.name}, [], 1);
cols.date  = reshape([entries.datenum], [], 1);
cols.bytes = reshape([entries.bytes], [], 1);

names = {};
for i = 1:n
    names = [names, reshape(fieldnames(entries(i).vals), 1, [])]; %#ok<AGROW>
end
names = [names, reshape(fieldnames(defaults), 1, [])];
if isempty(names), return, end
names = unique(names, 'stable');

for k = 1:numel(names)
    nm = names{k};
    raw = cell(n,1); has = false(n,1);
    for i = 1:n
        if isfield(entries(i).vals, nm)
            raw{i} = entries(i).vals.(nm); has(i) = true;
        elseif isfield(defaults, nm)
            raw{i} = defaults.(nm);        has(i) = true;
        end
    end
    present = raw(has);
    if all(cellfun(@(x) isnumeric(x) && isscalar(x), present))
        c = nan(n,1);
        c(has) = [present{:}];             % numeric column, NaN = missing
    elseif all(cellfun(@ischar, present))
        c = repmat({''}, n, 1);
        c(has) = present;                  % text column, '' = missing
    else
        c = raw;                           % vectors / mixed, [] = missing
    end
    cols.(nm) = c;
end
end
