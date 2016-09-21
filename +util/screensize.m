function varargout = screensize
%SCREENSIZE Return the size of current screen.
%   SZ = SCREENSIZE returns a vector [WIDTH, HEIGHT].
%   [W, H] = SCREENSIZE returns the width and height in W and H
%   respectively.

scnSize = get(0, 'ScreenSize');
scnSize = scnSize(3:4);

if nargout == 1
    varargout{1} = scnSize;
else
    varargout{1} = scnSize(1);
    varargout{2} = scnSize(2);
end

end

