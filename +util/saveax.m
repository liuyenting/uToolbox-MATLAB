function cdata = saveax(file, ax)
%SAVEAX Save an axes to TIFF file.
%   
%   SAVEAX(FILE) simply saves currently active axes to FILE, using TIFF
%   format.
%   SAVEAX(FILE, AX) saves the specified axes instead of the active one.
%
%   Note: Save operation is default to NEVER overwrite existing image. 

if nargin == 1
    ax = gca;
end

% extract the frame
frame = getframe(ax);
cdata = frame.cdata;

% save to file
tiff.imsave(cdata, file);

end

