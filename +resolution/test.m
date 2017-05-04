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
header = strrep(header, '"', '');
fclose(fid);

% data
data = csvread(filePath, 1, 0);
t = toc;

% result
fprintf('%d samples loaded, %.2fs elapsed\n', size(data, 1), t);
fprintf('%d columns in the dataset\n', length(header));

% identify target row
uncertInd = startsWith(header, 'uncertainty_xy');
uncertInd = find(uncertInd, 1, 'first');

fprintf('uncertainty column locates at %d\n', uncertInd);

uncert = data(:, uncertInd);
uncertAvg = mean(uncert);
uncertStd = std(uncert);

fprintf('mean = %.2f, s.d. = %.2f\n', uncertAvg, uncertStd);