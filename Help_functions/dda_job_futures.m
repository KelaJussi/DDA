function out = dda_job_futures(action, key, value)
% DDA_JOB_FUTURES  Keeps the parfeval futures of background DDA jobs.
% Used by dda_job_registry only.
%
%   dda_job_futures('set', key, future)
%   F = dda_job_futures('get', key)        [] if unknown
%   dda_job_futures('remove', key)
%
% Locked in memory with mlock so that "clear all" does not lose running
% background jobs. Because a locked function keeps its old code in memory
% even after the file is edited, this function is deliberately minimal and
% should never need changes; all other job logic lives in dda_job_registry,
% which is not locked.

persistent store
mlock
if isempty(store), store = containers.Map('KeyType','char','ValueType','any'); end
out = [];
key = char(key);
switch action
    case 'set',    store(key) = value;
    case 'get',    if isKey(store, key), out = store(key); end
    case 'remove', if isKey(store, key), remove(store, key); end
end
end
