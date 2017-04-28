function newcoord = offsetorigin(coord)
%OFFSETORIGIN Move the coordinates to starts from 0-based origin.
%
%   OFFSETORIGIN(COORD) returns the shifted coordinate row vector.

mincoord = min(coord);
n = size(coord, 1);
mincoord = repmat(mincoord, [n, 1]);
newcoord = coord - mincoord;

end
