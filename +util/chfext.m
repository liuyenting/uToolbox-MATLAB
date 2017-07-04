function newpath = chfext(oldpath, ext)
%CHFEXT Change file extension.
%   
%   NEWPATH = CHFEXT(OLDPATH, EXT) changes the file extension of OLDPATH to
%   EXT, if no file extension is detected, the extension is append
%   directly.

[fpath, fname, ~] = fileparts(oldpath);
newpath = fullfile(fpath, [fname, '.', ext]);

end

