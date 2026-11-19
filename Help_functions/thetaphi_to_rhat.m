function rhat = thetaphi_to_rhat(theta,phi)
% give theta and phi as matrices (from angexy_to_theta_phi()
% Returns 3D array of size (x,y,3), where 3rd dimension has the vector rhat
rhatx = sin(theta).*cos(phi);
rhaty = sin(theta).*sin(phi);
rhatz = cos(theta);
rhat = cat(3,rhatx,rhaty,rhatz);
end