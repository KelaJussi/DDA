function epsilon = epsilon_BK7(omega)
%omega in eV
    c = 2.99792458e8;
    hbar = 6.582119514e-16; 
    h = 2*pi*hbar;
    B1 = 1.03961212; B2 = 0.231792344; B3 = 1.01046945;
    A1 = 0.00600069867; A2 = 0.0200179144; A3 = 103.560653;
    C1 = h*c/sqrt(A1)/1e-6;C2 = h*c/sqrt(A2)/1e-6;C3 = h*c/sqrt(A3)/1e-6;
    D1 = B1*C1^2; D2 = B2*C2^2; D3 = B3*C3^2;

    epsilon = 1 + D1./(C1^2-omega.^2) + D2./(C2^2-omega.^2) + D3./(C3^2-omega.^2);
    
    %figure
    %plot(omega,sqrt(abs(epsilon)/2+real(epsilon)/2))
end