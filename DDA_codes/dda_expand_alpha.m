function Array = dda_expand_alpha(Array)
% DDA_EXPAND_ALPHA  Rebuild the full per-particle Array.alpha from the compact form.
% For the old solvers and builders (direct solver, legacy eigen solvers) that
% read Array.alpha(ii,:,:,kk) directly. Needs N times more memory than the
% compact form.
if ~isfield(Array, 'alpha_types'), return; end
A  = Array.alpha_types;
sz = size(A);
Array.alpha = reshape(A(Array.alpha_type, :), [numel(Array.alpha_type), sz(2:end)]);
Array = rmfield(Array, {'alpha_types', 'alpha_type'});
end
