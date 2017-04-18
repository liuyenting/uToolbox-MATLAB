function localizedemo
%LOCALDEMO Demonstrate the basic usage of class Localize.
%   LOCALDEMO will try to perform localization on presented image stack.

%% Load the stack.
imgStackPath = 'C:\Users\Lattice\Documents\MATLAB\data\local_seg_test\dense.tif';
tic;
[rawImgStack, rawStackSize] = tiff.imread(imgStackPath, true);
t = toc;
fprintf(' ** %f seconds to load the stack\n', t);
fprintf(' .. stack has the size of (x, y, z) = (%d, %d, %d)\n', ...
        rawStackSize(1), rawStackSize(2), rawStackSize(3));

end