function fwrite(outPath, type, arr, bdBox, comp)
%FWRITE Write as an Amira file
%   Detailed explanation goes here

%% Parse arguments
% enough arguments
if nargin < 3
    error('amira:fwrite', 'Not enough arguments');
elseif nargin > 5
    warning('amira:fwrite', 'Excess arguments are ignored.');
end

% verify input array dimension
nDim = ndims(arr);
if nDim > 3
    error('amira:fwrite', 'Array dimension should not exceed 3.');
end

% verify the bounding box
if nargin == 3
    bdBox = repmat([0, 1], [1, nDim]);
else
    if size(bdBox, 2)/2 ~= nDim
        error('amira:fwrite', 'Bounding box dimension mismatch, ignored.');
    end
end 

% verify compression type
if nargin < 5
    comp = 'Ascii'; %TODO: no compression, type string?
else
    switch(comp)
        case 'RLE'
            comp = 'HxByteRLE';
        case 'Zip'
            comp = 'HxZip';
        otherwise
            error('amira:fwrite', 'Unknown compressino scheme.');
    end
end
       
%% Create output file
if exist(outPath, 'file') == 2
    warning('amira:fwrite', 'Output file exists, overwrite.');
end

% open the file in binary mode
fid = fopen(outPath, 'w');
% register cleanup function upon exit.
hCleanup = onCleanup(@() fclose(fid));

% write header initializer
fprintf(fid, '# Avizo BINARY-LITTLE-ENDIAN 2.1\n\n');

%% Array info
arrSize = size(arr);

% write info to the file
fprintf(fid, 'define Lattice');
for s = arrSize
    fprintf(fid, ' %d', s);
end
fprintf(fid, '\n\n');

% write parameter declaration
fprintf(fid, 'Parameters {\n\tCoordType "uniform",\n');

% write bounding box
fprintf(fid, '\tBoundingBox');
for p = bdBox
    fprintf(fid, ' %f', p);
end
fprintf(fid, ',\n');

% write size definitions
fprintf(fid, '\tContent "%d', arrSize(1));
if length(arr) > 1
    for s = arrSize(2:end)
        fprintf(fid, 'x%d', s);
    end
end
% write data type string
typeStr = convType(arr);
fprintf(fid, ' %s', typeStr);

% write uniform definition
fprintf(fid, ', uniform coordinates"\n}\n\n');

% write Lattice type
fprintf(fid, 'Lattice { %s Labels } @1', typeStr);

% flatten the array
arr = arr(:);
% type cast to byte stream
arr = typecast(arr, 'uint8');
% compress the data if required
switch comp
    case { 'Ascii' }
        data = arr;
    case { 'HxByteRLE' }
        data = rleEncode(arr);
    case { 'HxZip' }
        error('amira:fwrite', 'Not implement yet.');
end
% count the length
len = length(data);

% write compression result
fprintf(fid, '(%s,%d) \n\n', comp, len); %TODO: write compression size
% write the data
fprintf(fid, '# Data section follows\n@1\n');
fwrite(fid, data);
fprintf(fid, '\n');

end

function typeStr = convType(obj)
%CONVTYPE Convert MATLAB data type to Amira data type.
%
%   TBA

% default type
typeStr = 'int8';

matlabType = {
    'uint8', 'single', 'int16', 'uint16', 'int32', 'double'
};
amiraType = {
    'byte', 'float', 'short', 'ushort', 'int', 'double'
};

for i = 1:length(matlabType)
    if isa(obj, matlabType{i})
        typeStr = amiraType{i};
        break;
    end
end

end

function out = rleEncode(in)
%RLEENCODE Encode data in RLE format.
%
%   TBA

% init to empty array
out = [];

nBytes = length(in);
i = 1;
while (i <= nBytes)
    if (i+1) > nBytes
        out = [out; 1; in(i)];
        i = i+1;
    elseif in(i) == in(i+1)
        j = 2;
        while (j < 127)
            if (i+j+1) > nBytes
                break
            elseif in(i) ~= in(i+j)
                break
            end
            j = j+1;
        end
        out = [out; j; in(i)];
        i = i+j;
    else
        j = 1;
        while (j < 127)
            if (i+j+1) > nBytes
                break
            elseif in(i+j) == in(i+j+1)
                break
            end
            j = j+1;
        end
        out = [out; bitor(j, 128, 'uint8'); in(i:i+j-1)];
        i = i+j;
    end
end

end
