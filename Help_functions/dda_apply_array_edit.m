function [Array, Input] = dda_apply_array_edit(Array, Input)
% DDA_APPLY_ARRAY_EDIT  Put the edited layout (DDA_array_editor) into the constructed array.
%
% Called by DDA_calc_looped right after Make_generic_Array. Input.array_edit
% (made by the array editor in the parameter GUI) holds:
%   r_i        N x 3 particle positions [m]
%   radius     N x 1 particle radii [m]; NaN = the radius from Input.particle.r
%   base_hash  fingerprint of the array the edits were made on
% The edits are applied only to that same array: if the inputs now construct a
% different array (other p, arrayN, lattice type ...), the run stops instead of
% putting positions meant for another array into it.
if ~isfield(Input, 'array_edit') || isempty(Input.array_edit) || ...
        ~isfield(Input.array_edit, 'enable') || ~Input.array_edit.enable
    return
end
E = Input.array_edit;
N = size(Array.r_i, 1);
if size(E.r_i, 1) ~= N || dda_checksum(Array.r_i) ~= string(E.base_hash)
    error('dda_apply_array_edit:mismatch', ['The array layout was edited for a different array ' ...
        '(%d particles then, %d now). The array inputs (p, arrayN, lattice type, shifts ...) changed ' ...
        'after the edit. Open the array editor again, or clear the layout edits.'], size(E.r_i,1), N);
end
if any(contains(string(Input.inversion_tool), "fft"))
    error('dda_apply_array_edit:fft', ['FFT-based solvers (%s) need the regular lattice; edited ' ...
        'layouts work with eig and the other solvers.'], Input.inversion_tool);
end
Array.r_i = cast(E.r_i, 'like', Array.r_i);
if isfield(Array, 'r_i_conv'), Array = rmfield(Array, 'r_i_conv'); end   % FFT grid, no longer valid
if isfield(E, 'radius') && any(~isnan(E.radius))
    r = E.radius(:);
    r(isnan(r)) = Input.particle.r;
    Array.radius = r;                          % per-particle radii, used by the polarizability
end
fprintf('Edited array layout applied: %d particle(s) moved, %d resized.\n', E.n_moved, E.n_resized);
end
