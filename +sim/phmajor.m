function B = phmajor(A, sz, no, np)
%PHMAJOR Convert input N-D array to orientation-phase major.
%   
%   B = PHMAJOR(A, SZ, NO, NP) reshapes A into orientation-phase majored B.
%   SZ is the size of a single stack, NO is the number of orientations
%   presented, NP is the phases presented.
%
%   Note
%   ----
%   All the orientations and phases are assumed to mix together in the Z
%   dimension, while orientation is the slowest dimension.
%   e.g. O1P1P2P3 O2P1P2P3 O3P1P2P3 (XYZ variations are not shown)

% flatten the original array
A = A(:);
% reduce total number of Zs
sz(3) = sz(3) / (no*np);

% reshape the dimension to consider orientations and phases
B = reshape(A, [sz(1:2), no, np, sz(3)]);

% re-order so that orientation and phases are the slowest
%   X Y O P Z  ->  O P Z X Y
%     X Y O P  ->  O P X Y
nd = length(sz);
if nd == 2
    order = [3, 4, 1, 2];
elseif nd == 3
    order = [3, 4, 5, 1, 2];
else
    error('sim:phmajor', 'Unsupported dimension format.');
end
B = permute(B, order);

end

