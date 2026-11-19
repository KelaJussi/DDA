function sec = dda_walltime_seconds(w)
% DDA_WALLTIME_SECONDS  Slurm wall time as Slurm applies it, in seconds.
%
%   sec = dda_walltime_seconds("0:30:00")     % 1800
%   sec = dda_walltime_seconds("0:00:10")     % 60: Slurm counts whole minutes
%
% Formats (as sbatch --time): m, m:s, h:m:s, d-h, d-h:m, d-h:m:s.
% Slurm keeps time limits in whole minutes and rounds seconds up, so a
% request of 10 s gives 1 min. NaN if unset or unreadable.
sec = NaN;
w = strtrim(string(w));
if ismissing(w) || strlength(w) == 0, return; end
days = 0;
if contains(w, "-")
    p = split(w, "-");
    days = str2double(p(1));
    f = str2double(split(p(2), ":")).';
    f = [f, zeros(1, 3 - numel(f))];                 % d-h, d-h:m, d-h:m:s
    s = days*86400 + f(1)*3600 + f(2)*60 + f(3);
else
    f = str2double(split(w, ":")).';
    switch numel(f)
        case 1, s = f(1)*60;                         % m
        case 2, s = f(1)*60 + f(2);                  % m:s
        case 3, s = f(1)*3600 + f(2)*60 + f(3);      % h:m:s
        otherwise, return
    end
end
if isnan(s), return; end
sec = 60*ceil(s/60);                                 % whole minutes, rounded up
end
