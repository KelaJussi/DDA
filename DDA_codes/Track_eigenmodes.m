function C = Track_eigenmodes(C, vectors, Input, Array)
% TRACK_EIGENMODES  Give every eigenvector a mode label (mx, my, polarization)
% and follow each labelled mode across the frequency grid.
%
% Replaces the mode-ordering code that was inlined in the old
% Calculate_{1D,2D,3D}_DDA_generic_eigenmodes files (fft_mode_search and
% mode_overlap_search) and produces the same fields JTK_plots reads.
%
% Inputs
%   C        from Calculate_DDA_eigenmodes: values, coefficients, exts_unordered
%   vectors  1 x Nk cell, eigenvectors V (Ndof x Ndof) of every frequency
%
% How a mode gets its label
%   Each eigenvector is projected onto separable basis functions evaluated at
%   the particle positions: sin((m+1)*pi*u) along open directions (the
%   standing-wave / DST basis; m = number of nodes) and exp(2i*pi*f*x/L) along
%   periodic directions. This is what the old padded-FFT peak search
%   approximated, but it also works for shifted or non-square lattices.
%   The fraction of the eigenvector's power that falls on each label is
%   its "basis score".
%
% How labels are assigned (Input.track.method)
%   "basis"   every frequency on its own: maximise the total basis score
%   "overlap" label the first point by basis score, then follow modes by the
%             overlap |<v_ref, v>|^2 with the nearest already-tracked point
%   "hybrid"  weighted sum of both costs (Input.track.overlap_weight)
%   Assignment is global (matchpairs), so two eigenvectors can never share a
%   label at the same frequency.
%
% On a complex-frequency grid the points are processed outward from the real
% axis, and each point is compared with the nearest point already tracked.
%
% Options, Input.track.* (defaults in brackets)
%   method          ["hybrid"]  "basis" | "overlap" | "hybrid"
%   overlap_weight  [0.5]       weight of the overlap cost in "hybrid"
%   basis           ["auto"]    "auto" (Fourier along periodic directions,
%                               sine otherwise) | "sine" | "fourier"
%   store_P_eigs    [true]      build C.P_eigs (Ndof x Nk per mode; large)
%   min_confidence  [0.5]       basis score below which a label counts as unclear
%   verbose         [true]
%
% Outputs added to C (shape = [nx ny], plus ncomp if > 1, plus nbranch if > 1)
%   C.vals            [shape Nk]  eigenvalue of each labelled mode
%   C.coefs_ordered   [shape Nk]  projection of Einc on each labelled mode
%   C.extsall         [Nk shape]  extinction of each labelled mode
%   C.P_eigs          cell(shape), each Ndof x Nk = v/lambda (single precision)
%   C.indices_ordered [shape Nk]  column of the unordered arrays for each label
%   C.track           method, label meaning, confidence, reference points
%
% DRAFT - written without MATLAB at hand.

opt = struct_defaults(getopt(Input,'track',struct()), struct( ...
    'method',"hybrid", 'overlap_weight',0.5, 'basis',"auto", ...
    'store_P_eigs',true, 'min_confidence',0.5, 'verbose',true));
method = lower(string(opt.method));
if ~ismember(method, ["basis","overlap","hybrid"])
    error('Track_eigenmodes:method','Input.track.method must be "basis", "overlap" or "hybrid".');
end

values = C.values;
[Ndof, Nk] = size(values);
N      = size(Array.r_i,1);
ncomp  = round(Ndof/N);
nx     = Input.arrayN(1);
ny     = 1; if numel(Input.arrayN) > 1, ny = Input.arrayN(2); end
ncell  = nx*ny;
nbranch = max(1, ceil(Ndof/(ncell*ncomp)));   % >1 e.g. for double_rectangular or stacked layers
nlab   = ncell*ncomp;                         % labels (mx,my,pol)
nslot  = nlab*nbranch;                        % labels x branches
slot_lab    = mod((1:nslot).'-1, nlab) + 1;   % label of each slot
slot_branch = floor(((1:nslot).'-1)/nlab) + 1;

%% Basis at the particle positions
Mb = dda_mode_basis(double(Array.r_i), nx, ny, opt.basis, Input);   % shared with the complex plots
lab_x = Mb.lab_x;  lab_y = Mb.lab_y;  basis_used = Mb.used;

%% Processing order and reference points
eV = Input.range.eV(:);
[~, order] = sortrows([abs(imag(eV)), real(eV)]);  % real axis first, then outward
ref = zeros(1,Nk);
for n = 2:Nk
    done = order(1:n-1);
    [~, i] = min(abs(eV(order(n)) - eV(done)));
    ref(order(n)) = done(i);
end

%% Assign labels point by point
idx  = NaN(nslot, Nk);          % slot -> column of the unordered arrays
conf = NaN(nslot, Nk);          % basis score of the assigned eigenvector
for n = 1:Nk
    k = order(n);
    V = vectors{k};
    score  = dda_mode_basis(Mb, V, ncomp);      % nlab x Ndof, columns sum to 1
    cost_b = 1 - score(slot_lab, :);            % nslot x Ndof

    r = ref(k);
    if r == 0 || method == "basis"
        cost = cost_b;
        use_branch_sort = true;
    else
        valid = ~isnan(idx(:,r));
        U = zeros(Ndof, nslot, 'like', V);
        U(:,valid) = vectors{r}(:, idx(valid,r));
        cost_o = 1 - abs(U'*V).^2;              % eigenvectors have unit norm
        cost_o(~valid,:) = cost_b(~valid,:);    % no reference for these slots
        if method == "overlap"
            cost = cost_o;
        else
            w = opt.overlap_weight;
            cost = w*cost_o + (1-w)*cost_b;
        end
        use_branch_sort = false;
    end

    M = matchpairs(max(cost,0), 1);             % costs are in [0,1], so every column is matched
    idx(M(:,1), k) = M(:,2);

    if use_branch_sort && nbranch > 1
        % Branches of one label have identical basis cost; order them by
        % Re(lambda) so the first labelled point is deterministic.
        for l = 1:nlab
            s  = find(slot_lab == l);
            js = idx(s,k); ok = ~isnan(js);
            [~, o] = sort(real(values(js(ok),k)));
            js_ok = js(ok); idx(s(ok),k) = js_ok(o);
        end
    end
    ok = ~isnan(idx(:,k));
    conf(ok,k) = score(sub2ind(size(score), slot_lab(ok), idx(ok,k)));
end

%% Build ordered outputs
shape = [nx ny];
if ncomp > 1 || nbranch > 1, shape(end+1) = ncomp; end
if nbranch > 1,              shape(end+1) = nbranch; end

vals  = complex(NaN(nslot,Nk));
coefs = complex(NaN(nslot,Nk));
exts  = NaN(nslot,Nk);
for k = 1:Nk
    ok = ~isnan(idx(:,k));
    j  = idx(ok,k);
    vals(ok,k)  = values(j,k);
    coefs(ok,k) = C.coefficients(j,k);
    exts(ok,k)  = C.exts_unordered(j,k);
end

C.vals            = reshape(vals,  [shape Nk]);
C.coefs_ordered   = reshape(coefs, [shape Nk]);
C.extsall         = permute(reshape(exts,[shape Nk]), [numel(shape)+1, 1:numel(shape)]);
C.indices_ordered = reshape(idx,   [shape Nk]);

if opt.store_P_eigs
    Pord = complex(NaN(Ndof, Nk, nslot, 'single'));   % single: half the file size, ample for plots
    for k = 1:Nk
        ok = find(~isnan(idx(:,k)));
        j  = idx(ok,k);
        Pord(:,k,ok) = reshape(vectors{k}(:,j)./values(j,k).', Ndof, 1, []);
    end
    C.P_eigs = reshape(num2cell(Pord,[1 2]), [shape 1]);
end

C.track = struct( ...
    'method', method, 'basis', basis_used, ...
    'labels_x', lab_x, 'labels_y', lab_y, ...           % what index 1..nx / 1..ny means
    'ncomp', ncomp, 'nbranch', nbranch, 'shape', shape, ...
    'confidence', reshape(conf,[shape Nk]), ...
    'n_unclear', sum(conf < opt.min_confidence, 1), ...
    'order', order(:).', 'reference', ref);

if opt.verbose
    fprintf(['Track_eigenmodes: %s, %s basis, %d modes x %d points; ' ...
             '%.1f%% of labels below confidence %.2f (worst point: %d unclear).\n'], ...
        method, basis_used, nslot, Nk, 100*mean(conf(:) < opt.min_confidence, 'omitnan'), ...
        opt.min_confidence, max(C.track.n_unclear));
end
end

%% ========================================================================
function s = struct_defaults(s, d)
f = fieldnames(d);
for n = 1:numel(f)
    if ~isfield(s,f{n}), s.(f{n}) = d.(f{n}); end
end
end

function v = getopt(s,name,default)
if isstruct(s) && isfield(s,name), v = s.(name); else, v = default; end
end
