function [C,Innz]=Make_checkerboard_array(nx,ny)
% makes a binary checkerboard matric of size n x n .
% For hexagonal arrays, the size of the array of scatterers (using fft
% speed-up for convolution) has to be of size 2n x n to keep the number
% of active scattering particles the same as with rectangular arrays.
%
%  x . x . x . x .   , where x = location of the active scatterer
%  . x . x . x . x   ,       . = location of the fft dipole site
%  x . x . x . x .
%  . x . x . x . x
%
% Real array size is thus 2*nx x ny  

%% builds the checkerboard array C, using eye(2) as a building block:
% C_tmp=repmat(eye(2),[nx ceil(ny/2)]);
num_zero = 1e-3;    % numerical zero
% C_tmp=repmat([1 num_zero ; num_zero 1],[nx ceil(ny/2)]);
C_tmp=repmat([1 num_zero ; num_zero 1],[nx ny]);

%% remove the additional line if ny is odd, since the eye(2) always makes a even array:
if mod(ny,2) == 0 
    C=C_tmp; 
else  
    C=C_tmp(1:(end-1),1:(end-1));
end

%% find the non-zero indices to save memory:
Innz=find(C);

end