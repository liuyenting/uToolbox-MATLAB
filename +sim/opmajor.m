function varargout = opmajor(I, no, np)
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

%% parameters
sz = size(I);
% size of a single plane
imSz = sz(1:2);
% number of slices (Z)
nz = sz(3) / (no*np);

%% process
% reshape the dimension to consider orientations and phases
I = reshape(I, [imSz, np, no, nz]);

% re-order the stack so that orientations and phases are the slowest
%   X Y P O Z  ->  X Y | Z P O
%     X Y O P  ->  X Y | 1 P O
% maximum order hard-coded to 5 indicates that maximum supported dimension
% is 3-D
order = [1, 2, 5, 3, 4];
I = permute(I, order);

%% output
if nargout == 1
    varargout{1} = I;
elseif nargout == 2
    varargout{2} = size(I);
end

end

