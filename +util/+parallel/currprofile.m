function profile = currprofile
%CURRPROFILE Get current running profile.

poolobj = gcp('nocreate');
if isempty(poolobj)
    profile = 'none';
else
    profile = poolobj.Cluster.Profile;
end

end

