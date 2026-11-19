function Input = dda_update_derived_inputs(Input, before, use_eV, wl_to_eV)
% DDA_UPDATE_DERIVED_INPUTS  Recompute what Make_predefined_Inputs derives from
% options the parameter GUI (or saved defaults) can change, using the same
% formulas, but only when those options changed. Used by Make_predefined_Inputs
% (after saved defaults and the GUI) and by the array editor (to build the array
% from the inputs as they are in the GUI).
%
%   before    the Input before the changes
%   use_eV    false if Make_predefined_Inputs defined the range in wavelength
%   wl_to_eV  1239.84187e-9 (wavelength in m times energy in eV)
if nargin < 3, use_eV = true; end
if nargin < 4, wl_to_eV = 1239.84187e-9; end
if ~isequal(Input.range.real_parts, before.range.real_parts) || ...
        ~isequal(Input.range.imag_parts, before.range.imag_parts)
    if ~use_eV
        disp('The GUI range is in eV: using range.real_parts/imag_parts instead of wl.')
    end
    Input.range.eV  = reshape(Input.range.real_parts + Input.range.imag_parts.'*1i, 1, ...
                              length(Input.range.real_parts)*length(Input.range.imag_parts));
    Input.range.wl  = wl_to_eV./Input.range.eV;
    Input.range.k   = 2*pi./Input.range.wl;
    Input.range.k_N = length(Input.range.k);
end
if Input.array_type == "double_rectangular"
    Input.nmbofmodes = Input.arrayN(1)*Input.arrayN(2)*2;
else
    Input.nmbofmodes = Input.arrayN(1)*Input.arrayN(2);
end
if ~isequal(Input.p, before.p) || ~isequal(Input.arrayN, before.arrayN)
    % follow p and arrayN, unless the value itself was edited in the GUI
    if isequal(Input.second_array_offset, before.second_array_offset)
        Input.second_array_offset = [Input.p(1)/4, Input.p(2)/4, 0];
    end
    if isequal(Input.Gaussian_sigma, before.Gaussian_sigma)
        Input.Gaussian_sigma = Input.arrayN(1)*Input.p(1)/4*1.5;
    end
end
if Input.array_type == "PC_mirror_reflection_test"
    Input.preset_field = 'PC_mirror_reflection_test';     % as in Make_predefined_Inputs
end
end
