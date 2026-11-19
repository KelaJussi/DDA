function h = dda_checksum(varargin)
% DDA_CHECKSUM  Short hexadecimal fingerprint of numeric arrays.
%   h = dda_checksum(A, B, ...)   e.g. dda_checksum(r_i) or dda_checksum(r_i, radius)
% Values are rounded to 1e-12 (picometres for positions in metres), so the same
% layout always gives the same fingerprint. Used for array-layout edits: in the
% result file name, and to check that edits still belong to the array the
% inputs construct. NaN values count as "not set".
s1 = 0; s2 = 0; n0 = 0;
for a = varargin
    x = double(a{1}(:));
    x(isnan(x)) = -1;
    q = round(x*1e12);
    k = (1:numel(q)).' + n0;
    s1 = mod(s1 + sum(mod(q .* mod(k*2654435761, 1000003), 4294967291)), 4294967291);
    s2 = mod(s2 + sum(mod(q + 7919*k, 4294967279)), 4294967279);
    n0 = n0 + numel(q);
end
h = sprintf('%08x%04x', mod(s1, 2^32), mod(s2, 2^16));
end
