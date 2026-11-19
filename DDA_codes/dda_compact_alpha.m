function Array = dda_compact_alpha(Array)
% DDA_COMPACT_ALPHA  Store one polarizability per particle type instead of per particle.
%
% Replaces Array.alpha (N x 3 x 3 x Nk, or N x Nk for cyl_TM) by
% Array.alpha_types (M x ...) and Array.alpha_type (N x 1), where M is the
% number of distinct particles: 1 for identical particles, a few for wall or
% Gaussian lattices. Read it with dda_alpha; dda_expand_alpha undoes this.
if ~isfield(Array, 'alpha') || isfield(Array, 'alpha_types'), return; end
A  = Array.alpha;
sz = size(A);  N = sz(1);
R  = reshape(A, N, []);
if all(all(R == R(1,:)))                     % identical particles (the usual case)
    idx = ones(N, 1);  first = 1;
else
    % group equal rows by a random projection (equal rows give equal keys),
    % then verify exactly; fall back to a direct comparison if keys collide
    key = [real(R) imag(R)] * randn(2*size(R,2), 1);
    [~, first, idx] = unique(key, 'stable');
    if ~isequal(R, R(first(idx), :))
        first = zeros(0,1); idx = zeros(N,1);
        for n = 1:N
            m = find(all(R(first,:) == R(n,:), 2), 1);
            if isempty(m), first(end+1,1) = n; m = numel(first); end %#ok<AGROW>
            idx(n) = m;
        end
    end
end
Array.alpha_types = reshape(R(first, :), [numel(first), sz(2:end)]);
Array.alpha_type  = idx(:);
Array = rmfield(Array, 'alpha');
end
