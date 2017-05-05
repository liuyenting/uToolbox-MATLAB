clear all; close all; %#ok<CLALL>

fprintf('\n -- loading the data --\n');

filePath = fullfile(userpath, '0421Neuron1_postprocessed.csv');

tic;
% header
fid = fopen(filePath);
header = fgetl(fid);
% split by comma
header = strsplit(header, ',');
% remove quotes
header = strrep(header, '"', '')
fclose(fid);

% data
%data = csvread(filePath, 1, 0);
t = toc;

% result
%fprintf('%d samples loaded, %.2fs elapsed\n', size(data, 1), t);
fprintf('%d columns in the dataset\n', length(header));

% identify target row
query = {'uncertainty_xy', 'x', 'y'};
uncertInd = findcol(header, query)

fprintf('uncertainty column locates at %d\n', uncertInd);
