function newcoord = offsetorigin(coord)
%OFFSETORIGIN Move the coordinates to starts 0-based origin.
%
%   OFFSETORIGIN(COORD) returns the shifted coordinate row vector.

mincoord = min(coord);

% pad the subtraction factor
n = size(coord, 1);
mincoord = repmat(mincoord, [n, 1]);

% shift the coordinates
newcoord = coord - mincoord;

end
