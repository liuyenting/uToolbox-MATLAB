function [B, nsz] = opmajor(A, osz, no, np)
%OPMAJOR Convert input N-D array to orientation-phase major.
%   
%   [B, NSZ] = OPMAJOR(A, OSZ, NO, NP) reshapes A into orientation-phase 
%   majored B. OSZ is the size of a single stack, NO is the number of 
%   orientations presented, NP is the phases presented. Updated volume
%   size, which takes consider NO and NP is saved as NSZ.
%
%   Note
%   ----
%   All the orientations and phases are assumed to mix together in the Z
%   dimension, while orientation is the slowest dimension.
%   e.g. O1P1P2P3 O2P1P2P3 O3P1P2P3 (XYZ variations are not shown)

% flatten the original array
A = A(:);
% reduce total number of Zs
osz(3) = osz(3) / (no*np);

% reshape the dimension to consider orientations and phases
B = reshape(A, [osz(1:2), no, np, osz(3)]);

% re-order so that orientation and phases are the slowest
%   X Y O P Z  ->  X Y | Z P O
%     X Y O P  ->  X Y | 1 P O
% maximum order hard-coded to 5 indicates that maximum supported dimension
% is XYZ
order = 1:5;
order = [order(1:2), order(end:-1:3)];
B = permute(B, order);

% save the new volume size
nsz = osz;

end

