function report = compare_mode_tracking(varargin)
% COMPARE_MODE_TRACKING  Compare legacy and new eigenmode labels of the same run.
%
%   report = compare_mode_tracking(Input)
%       Uses the two real-frequency result files of this Input: the normal
%       one and the one saved with Input.use_legacy_eig_solvers = true
%       (" legacy" before .mat), both in dda_results_folder("real").
%   report = compare_mode_tracking(file_new, file_legacy)
%   report = compare_mode_tracking(..., 'Overlap', 0.999, 'Plot', true, 'MaxRows', 20)
%
% How labels are compared
%   Neither code stores the ordered eigenvalues in the same form, but both
%   store C.P_eigs{label}(:,k) = v/lambda. For each frequency, every legacy
%   column is matched to the new column it overlaps most with:
%       overlap = |<a,b>| / (|a| |b|)      (1 = same eigenvector, any phase)
%   and |a|/|b| = |lambda_new|/|lambda_legacy| (1 = same eigenvalue).
%   The legacy label then "agrees" if the new code put that same
%   eigenvector on the same (mx,my) label.
%
% Outcome per legacy label and frequency
%   agree       same eigenvector, same (mx,my) label
%   relabelled  same eigenvector, different label in the new code
%   no_match    no new eigenvector overlaps >= 'Overlap' (the two runs
%               solved different matrices, or the new file lacks that column)
%   empty       the legacy code left this label unassigned (NaN)
%   duplicate   (flag) the legacy code put this eigenvector on two labels
%
% Polarization: the legacy 3D solver stores P_eigs without a polarization
% index, so only (mx,my) is compared there. The legacy 2D solver stores
% {mx,my,pol}; the polarization match is then reported separately.
% Periodic directions: the new labels are signed spatial frequencies, the
% legacy ones are magnitudes, so |label| is compared there.
%
% Output (also printed): report.summary, report.counts, report.per_label,
% report.disagreements (a table in MATLAB), report.status [legacy labels x Nk].

%% Arguments
[files, opt] = parse_args(varargin{:});
L = load_side(files{2}, 'legacy');
N = load_side(files{1}, 'new');
if L.Nk ~= N.Nk || L.Ndof ~= N.Ndof
    error('compare_mode_tracking:size', ...
        'The files differ in size: legacy %d dofs x %d frequencies, new %d x %d.', ...
        L.Ndof, L.Nk, N.Ndof, N.Nk);
end
lam = N.Input.range.wl(:).'*1e9;

%% New labels, in legacy-comparable form
[nx_new, ny_new, pol_new] = new_labels(N);       % one entry per new slot
[mx_leg, my_leg, pol_leg] = legacy_labels(L);    % one entry per legacy slot
has_pol = ~isempty(pol_leg);

%% Match eigenvectors frequency by frequency
codes = struct('agree',1, 'relabelled',2, 'no_match',3, 'empty',4);
nL = numel(L.Pe);  Nk = N.Nk;
status    = zeros(nL, Nk);
best_ov   = NaN(nL, Nk);  ratio = NaN(nL, Nk);  match = NaN(nL, Nk);
dup       = false(nL, Nk);
new_only  = zeros(1, Nk);
for k = 1:Nk
    A = column_block(L.Pe, k);   B = column_block(N.Pe, k);
    va = all(isfinite(A),1) & any(A,1);
    vb = all(isfinite(B),1) & any(B,1);
    status(~va, k) = codes.empty;
    if ~any(va), continue; end
    if ~any(vb), status(va, k) = codes.no_match; continue; end
    An = A(:,va)./vecnorm(A(:,va));  Bn = B(:,vb)./vecnorm(B(:,vb));
    [ov, j] = max(abs(An'*Bn), [], 2);
    ib = find(vb);  j = ib(j);
    ia = find(va);
    best_ov(ia,k) = ov;
    match(ia,k)   = j;
    ratio(ia,k)   = vecnorm(A(:,ia)) ./ vecnorm(B(:,j));
    ok = ov >= opt.Overlap;
    same = ok & mx_leg(ia) == nx_new(j) & my_leg(ia) == ny_new(j);
    status(ia(same), k)        = codes.agree;
    status(ia(ok & ~same), k)  = codes.relabelled;
    status(ia(~ok), k)         = codes.no_match;
    % the same new eigenvector matched by two legacy labels
    jj = j(ok);  [~, ~, g] = unique(jj);  cnt = accumarray(g(:), 1);
    dup(ia(ok), k) = cnt(g) > 1;
    new_only(k) = numel(setdiff(find(vb), unique(jj)));
end

%% Summary
assigned = status ~= codes.empty;
n_assigned = nnz(assigned);
c.agree      = nnz(status == codes.agree);
c.relabelled = nnz(status == codes.relabelled);
c.no_match   = nnz(status == codes.no_match);
c.empty      = nnz(status == codes.empty);
c.duplicate  = nnz(dup);
c.new_only   = sum(new_only);
report.counts = c;

ok_ov = best_ov(status == codes.agree | status == codes.relabelled);
report.summary = sprintf(['%d legacy label-frequency pairs assigned: %.1f%% agree, %.1f%% relabelled, ' ...
    '%.1f%% without a matching eigenvector. Legacy duplicates: %d. Labels left empty by legacy: %d; ' ...
    'eigenvectors labelled only by the new code: %d.'], n_assigned, 100*c.agree/max(n_assigned,1), ...
    100*c.relabelled/max(n_assigned,1), 100*c.no_match/max(n_assigned,1), c.duplicate, c.empty, c.new_only);

fprintf('\n== Mode tracking: legacy vs new ==\n  legacy: %s\n  new:    %s\n\n%s\n', files{2}, files{1}, report.summary);
if ~isempty(ok_ov)
    r = ratio(status == codes.agree | status == codes.relabelled);
    fprintf('Matched eigenvectors: median overlap %.6f, |lambda| ratio within %.1e of 1.\n', ...
        median(ok_ov), max(abs(r - 1)));
end
if ~has_pol && ~isempty(pol_new)
    fprintf(['The legacy file stores one polarization per (mx,my) (legacy 3D), so most of the\n' ...
             'eigenvectors "labelled only by the new code" are simply the other polarizations.\n']);
end
if c.no_match > 0.2*max(n_assigned,1)
    fprintf(['Note: many legacy eigenvectors have no counterpart. The two runs probably solved\n' ...
             'different matrices (e.g. periodic boundaries or mirrors the legacy builders lack),\n' ...
             'so their labels cannot be compared one to one.\n']);
end

%% Disagreements, per label and in detail
lab = arrayfun(@(i) label_str(mx_leg(i), my_leg(i), pol_leg, i), 1:nL, 'UniformOutput', false);
per = struct('label', lab, ...
    'agree',      num2cell(sum(status == codes.agree, 2)).', ...
    'relabelled', num2cell(sum(status == codes.relabelled, 2)).', ...
    'no_match',   num2cell(sum(status == codes.no_match, 2)).');
[~, order] = sort([per.relabelled] + [per.no_match], 'descend');
report.per_label = per(order);

[li, ki] = find(status == codes.relabelled | status == codes.no_match);
rows = struct('wavelength_nm', {}, 'legacy_label', {}, 'new_label', {}, 'overlap', {}, ...
              'new_confidence', {}, 'legacy_duplicate', {});
for r = 1:numel(li)
    i = li(r); k = ki(r); j = match(i,k);
    if isnan(j), nl = '-'; cf = NaN;
    else, nl = label_str(nx_new(j), ny_new(j), pol_new, j); cf = N.conf(j,k);
    end
    rows(end+1) = struct('wavelength_nm', lam(k), 'legacy_label', lab{i}, 'new_label', nl, ... %#ok<AGROW>
        'overlap', best_ov(i,k), 'new_confidence', cf, 'legacy_duplicate', dup(i,k));
end
if exist('struct2table','file') && ~isempty(rows), report.disagreements = struct2table(rows);
else, report.disagreements = rows; end

if has_pol
    same_pol = 0; both = 0;
    for k = 1:Nk
        i = find(status(:,k) == codes.agree);
        both = both + numel(i);
        same_pol = same_pol + nnz(pol_leg(i) == pol_new(match(i,k)));
    end
    fprintf('Polarization index also agrees for %.1f%% of the agreeing labels.\n', 100*same_pol/max(both,1));
end

if ~isempty(rows)
    fprintf('\nLabels with most disagreements:\n');
    for q = 1:min(8, numel(report.per_label))
        p = report.per_label(q);
        if p.relabelled + p.no_match == 0, break; end
        fprintf('  %-12s relabelled %4d, no match %4d, agree %4d\n', p.label, p.relabelled, p.no_match, p.agree);
    end
    fprintf('\nFirst %d disagreements (new_confidence = basis score of the new label, 0..1):\n', ...
        min(opt.MaxRows, numel(rows)));
    fprintf('  %10s  %-12s %-12s %8s %8s %s\n', 'wl (nm)', 'legacy', 'new', 'overlap', 'conf', 'dup');
    for r = 1:min(opt.MaxRows, numel(rows))
        w = rows(r);
        fprintf('  %10.2f  %-12s %-12s %8.4f %8.3f %d\n', w.wavelength_nm, w.legacy_label, ...
            w.new_label, w.overlap, w.new_confidence, w.legacy_duplicate);
    end
    lowconf = [rows.new_confidence];
    fprintf('\nMedian new confidence: %.3f at disagreements vs %.3f at agreements.\n', ...
        median(lowconf, 'omitnan'), median(conf_at(N, match, status == codes.agree), 'omitnan'));
else
    fprintf('No disagreements.\n');
end
report.status = status;  report.codes = codes;  report.labels = lab;
report.wavelength_nm = lam;  report.overlap = best_ov;  report.files = files;

%% Plots
if opt.Plot, plot_report(report, L, N, mx_leg, my_leg, pol_leg, codes); end
end

%% ========================================================================
function [files, opt] = parse_args(varargin)
opt = struct('Overlap', 0.999, 'Plot', true, 'MaxRows', 20);
if nargin >= 1 && isstruct(varargin{1})
    Input = varargin{1};
    Input.use_legacy_eig_solvers = false;
    name  = construct_filename(Input, "real");
    folder = dda_results_folder("real");
    files = {fullfile(folder, name), fullfile(folder, replace(name, ".mat", " legacy.mat"))};
    rest  = varargin(2:end);
elseif nargin >= 2
    files = {char(varargin{1}), char(varargin{2})};
    rest  = varargin(3:end);
else
    error('Give an Input struct, or the new and the legacy result files.');
end
for n = 1:2:numel(rest), opt.(rest{n}) = rest{n+1}; end
for n = 1:2
    if ~isfile(files{n}) && isfile(fullfile(dda_results_folder("real"), files{n}))
        files{n} = fullfile(dda_results_folder("real"), files{n});
    end
    if ~isfile(files{n}), error('compare_mode_tracking:file', 'File not found: %s', files{n}); end
end
end

function S = load_side(file, which)
D = load(file, 'C', 'Input');
if ~isfield(D.C, 'P_eigs')
    error('compare_mode_tracking:noPeigs', ...
        ['The %s file has no C.P_eigs. Rerun with Input.save.P_eigs = true ' ...
         '(and a larger Input.save.max_file_GB if it was dropped for size).'], which);
end
S.Input = D.Input;
S.C     = D.C;
S.Pe    = D.C.P_eigs(:);
S.shape = size(D.C.P_eigs);
first   = find(~cellfun(@isempty, S.Pe), 1);
[S.Ndof, S.Nk] = size(S.Pe{first});
if isfield(D.C, 'track')
    S.conf = reshape(D.C.track.confidence, [], S.Nk);
else
    S.conf = NaN(numel(S.Pe), S.Nk);
end
end

function [mx, my, pol] = new_labels(N)
if isfield(N.C, 'track')
    shape = N.C.track.shape; lx = N.C.track.labels_x; ly = N.C.track.labels_y;
else
    shape = N.shape; lx = 0:shape(1)-1; ly = 0:shape(2)-1;
end
sub = cell(1, numel(shape)); [sub{:}] = ind2sub([shape 1], (1:prod(shape)).');
mx = abs(lx(sub{1})); my = abs(ly(sub{2}));      % |label|: Fourier labels are signed
mx = mx(:); my = my(:);
pol = [];
if numel(shape) >= 3 && shape(3) > 1, pol = sub{3}; end
end

function [mx, my, pol] = legacy_labels(L)
sub = cell(1, max(2, numel(L.shape))); [sub{:}] = ind2sub([L.shape 1], (1:numel(L.Pe)).');
mx = sub{1} - 1;  my = sub{2} - 1;               % legacy labels are 0-based node counts
pol = [];
if numel(L.shape) >= 3 && L.shape(3) > 1, pol = sub{3}; end
end

function M = column_block(Pe, k)
M = complex(NaN(size(Pe{find(~cellfun(@isempty,Pe),1)},1), numel(Pe)));
for s = 1:numel(Pe)
    if ~isempty(Pe{s}), M(:,s) = Pe{s}(:,k); end
end
end

function s = label_str(mx, my, pol, i)
s = sprintf('(%d,%d)', mx, my);
if ~isempty(pol), names = 'xyz'; s = sprintf('%s %s', s, names(min(pol(i),3))); end
end

function v = conf_at(N, match, mask)
[i, k] = find(mask);
v = NaN(numel(i), 1);
for r = 1:numel(i), v(r) = N.conf(match(i(r),k(r)), k(r)); end
end

function plot_report(R, L, N, mx_leg, my_leg, pol_leg, codes)
dda_figure([100 100 1200 500], 'Name','mode tracking: legacy vs new');
tl = tiledlayout(1, 2, 'TileSpacing','compact');

nexttile(tl)
imagesc(R.wavelength_nm, 1:numel(R.labels), R.status, [1 4]);
colormap(gca, [0.30 0.70 0.40; 0.90 0.55 0.15; 0.80 0.20 0.20; 0.92 0.92 0.92]);
cb = colorbar; cb.Ticks = [1.375 2.125 2.875 3.625];
cb.TickLabels = {'agree','relabelled','no match','empty'};
nl = numel(R.labels); step = max(1, ceil(nl/25));
yticks(1:step:nl); yticklabels(R.labels(1:step:nl));
xlabel('Wavelength (nm)'); ylabel('Legacy label'); title('Legacy label vs new label')

nexttile(tl); hold on
[~, worst] = sort(sum(R.status == codes.relabelled | R.status == codes.no_match, 2), 'descend');
cols = lines(3); leg = {};
for q = 1:min(3, nl)
    i = worst(q);
    if all(R.status(i,:) == codes.agree | R.status(i,:) == codes.empty), break; end
    al = arrayfun(@(k) norm(L.Pe{i}(:,k)), 1:numel(R.wavelength_nm));
    [nxn, nyn, poln] = new_labels(N);
    s = find(nxn == mx_leg(i) & nyn == my_leg(i));          % same (mx,my) in the new labels
    if ~isempty(pol_leg) && ~isempty(poln), s = s(poln(s) == pol_leg(i)); end
    an = NaN(numel(s), numel(R.wavelength_nm));
    for t = 1:numel(s), an(t,:) = arrayfun(@(k) norm(N.Pe{s(t)}(:,k)), 1:numel(R.wavelength_nm)); end
    plot(R.wavelength_nm, al, '-', 'Color', cols(q,:), 'LineWidth', 2)
    plot(R.wavelength_nm, max(an, [], 1), '--', 'Color', cols(q,:), 'LineWidth', 1.5)   % largest over pol/branch
    leg(end+1:end+2) = {[R.labels{i} ' legacy'], [R.labels{i} ' new']}; %#ok<AGROW>
end
hold off; box on; set(gca, 'YScale', 'log')
if ~isempty(leg), legend(leg, 'Location', 'best'); end
xlabel('Wavelength (nm)'); ylabel('1/|\lambda|  (= |P_{eig}|)')
title('Labels with most disagreements')
end
