function tf = dda_alpha_is_SI(Input, Array)
% DDA_ALPHA_IS_SI  True when Array.alpha is an SI polarizability (C m^2/V).
%
% The SI convention of the solver (Make_int_matrix_A_SI and the Mie path):
%   alpha_SI = 4*pi*eps0*n_h^2 * alpha_vol,   p = alpha_SI * E
%   A = diag(1/alpha_SI) - k^2/(eps0*n_h^2) * G,  k = n_h*k0
% The old Gaussian-form builders (Make_*_int_matrix_A_sub*, Make_Afft_*) use
% exp(ikR)/R without these factors, so they must divide their interaction
% terms by 4*pi*eps0*n_h^2 whenever the polarizabilities are SI, and the
% extinction formulas need the same factor.
%
% Replaces the sample-type name checks ("Au_Mie" || "glass_Mie" || ...) that
% were copied, with different lists, into the builders and solvers.
% Polarizability functions that produce SI values mark it with
% Array.alpha_units = "SI"; the Mie sample types are SI by construction.

tf = any(string(Input.sample.type) == ["Au_Mie","glass_Mie","Au_Drude_Mie"]);
if nargin > 1 && isstruct(Array) && isfield(Array,'alpha_units')
    tf = string(Array.alpha_units) == "SI";
end
end
