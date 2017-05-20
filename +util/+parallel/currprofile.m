function profile = currprofile
%CURRPROFILE Get current running profile.
%
%   PROFILE = CURRPROFILE() returns the current running cluster profile,
%   'none' is returned if no parallel pool is running.

poolobj = gcp('nocreate');
if isempty(poolobj)
    profile = 'none';
else
    profile = poolobj.Cluster.Profile;
end

end

