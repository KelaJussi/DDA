function ratio = Calculate_mode_absorption_ratio(mode,k0,A_offdiag,A_offdiag_real, n_host)
% k0 is real(k0), as the first order approximation
    if nargin == 4
        n_host = 1.5;
    end
    c0 = 2.99979e8;
    eps0 = 8.854e-12; 
    omega_r = k0 * c0; 
    p = reshape(mode, [],1);

% P_ext, unchanged from before — plain transpose .'
Eloc = - double(A_offdiag) * p;
P_ext_tot = (omega_r/2)*imag(p.'*conj(Eloc));

% P_rad, now including interference — conjugate transpose '
P_rad_independent = sum(n_host*omega_r^4*abs(p).^2)/(12*pi*eps0*c0^3);
P_rad_interference = -(omega_r/2)*imag(p'*double(A_offdiag_real)*p);
P_rad_tot = P_rad_independent + P_rad_interference;

P_abs_tot = P_ext_tot - P_rad_tot;
ratio = P_abs_tot / P_ext_tot;




%    Pabs = omega/2 * sum( imag( mode .* conj( alpha.^(-1) .* mode) ) - 2/3 * k^3 * abs(mode).^2,"all");
 %   Pext = omega/2 * sum( imag( mode .* conj( alpha.^(-1) .* mode) ), "all");
  %  ratio = Pabs / Pext;

end