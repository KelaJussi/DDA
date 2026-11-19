function [theta,phi] = anglexy_to_theta_phi(x,y)
%converts cartesian to polar
% Have x as row vector [1,x], y as column vector [y,1]. 
%Outputs [y,x] matrices
% use radians
    theta = sqrt(x.^2+y.^2);
    phi = atan2(y,x);
end
