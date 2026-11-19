function ok = dda_cluster_reconnect(profiles)
% DDA_CLUSTER_RECONNECT  Renew MATLAB's SSH connection to the cluster without restarting.
%
%   dda_cluster_reconnect                       all remote cluster profiles
%   dda_cluster_reconnect("Roihu R2026b")       one profile
%   ok = dda_cluster_reconnect(...)             true if every connection works again
%
% When the SSH certificate expires (CSC: after 24 h), the connection MATLAB's
% cluster plugin opened earlier dies, but MATLAB keeps reusing it. Submitting
% then fails with "Failed to start mirror ... Unable to send channel-open
% request", even though a new certificate works for a fresh connection.
% Restarting MATLAB discards the dead connection; this function replaces it:
%   1. clears the places that hold on to the old connection: the plugin's
%      getRemoteConnection (if it caches in a persistent variable) and the
%      cluster objects cached by dda_job_registry (background jobs are kept,
%      they live in dda_job_futures);
%   2. reconnects the connection object the plugin uses (asked from the plugin's
%      own getRemoteConnection; otherwise MATLAB's cached object for the host,
%      user and remote job storage); reconnect recreates its SSH connections;
%   3. tests it with a remote command;
%   4. deletes jobs left in the "pending" state by failed submissions (they never
%      reached Slurm).
% First renew the certificate (MyCSC or CSC's certificate tool), so that it is
% where the profile looks for it (SSH agent / Pageant, or next to the key file).
% submit_dda_job calls this automatically once when a submission fails this way.

if nargin < 1 || isempty(profiles)
    profiles = string(parallel.clusterProfiles());
end
profiles = string(profiles);

% 1. forget cached connection holders
clear('getRemoteConnection');            % plugin function, if it caches in a persistent variable
clear('dda_job_registry');               % its cached cluster objects (not locked; futures are kept)

ok = true;
n_remote = 0;
for p = profiles(:).'
    try
        cl = parcluster(char(p));
    catch
        continue
    end
    if ~isprop(cl, 'AdditionalProperties'), continue; end
    ap = cl.AdditionalProperties;
    host = prop(ap, 'ClusterHost');  user = prop(ap, 'Username');
    store = prop(ap, 'RemoteJobStorageLocation');
    if strlength(host) == 0, continue; end                 % not a remote-submission profile
    n_remote = n_remote + 1;
    fprintf('Profile "%s" (%s@%s):\n', p, user, host);

    % 2. reconnect the connection object the plugin uses: ask the plugin itself
    %    (its getRemoteConnection), otherwise MATLAB's cache of connections
    try
        r = plugin_connection(cl);
        if isempty(r) && strlength(store) > 0
            r = parallel.cluster.RemoteClusterAccess.getConnectedAccessWithMirror( ...
                    char(host), char(store), char(user));
        elseif isempty(r)
            r = parallel.cluster.RemoteClusterAccess.getConnectedAccess(char(host), char(user));
        end
        reconnect(r);
        % 3. test it
        [status, out] = runCommand(r, 'echo connection-ok');
        if status == 0 && contains(string(out), "connection-ok")
            fprintf('  connection renewed and working.\n');
        else
            ok = false;
            fprintf('  reconnected, but the test command failed (status %d): %s\n', status, strtrim(out));
        end
    catch err
        ok = false;
        fprintf('  could not reconnect: %s\n', err.message);
        fprintf(['  Check that the new certificate is in place, or connect once with your own\n' ...
                 '  RemoteClusterAccess object; if it still fails, restarting MATLAB is the fallback.\n']);
        continue
    end

    % 4. clean up jobs that failed during submission
    try
        stuck = findJob(cl, 'State', 'pending');
        if ~isempty(stuck)
            ids = [stuck.ID];
            delete(stuck);
            fprintf('  deleted %d job(s) left from failed submissions: %s\n', numel(ids), mat2str(ids));
        end
    catch err
        fprintf('  (could not check for jobs left from failed submissions: %s)\n', err.message);
    end
end
if n_remote == 0
    fprintf('No remote cluster profiles found.\n');
    ok = false;
end
if nargout == 0, clear ok; end
end

function r = plugin_connection(cl)
% The connection object of the cluster plugin (getRemoteConnection in the
% plugin folder), or [] if the plugin has no such function
r = [];
pdir = '';
try, pdir = char(cl.PluginScriptsLocation); catch, end
if isempty(pdir) || ~isfile(fullfile(pdir, 'getRemoteConnection.m')), return; end
old = addpath(pdir);
try
    r = getRemoteConnection(cl);
catch
    r = [];
end
path(old);
end

function v = prop(ap, name)
% one AdditionalProperties value as a string ("" if missing)
v = "";
try
    if isprop(ap, name) || (isstruct(ap) && isfield(ap, name)), v = string(ap.(name)); end
catch
    try, v = string(ap.(name)); catch, end
end
if ismissing(v), v = ""; end
end
