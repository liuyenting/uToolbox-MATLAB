%close all;
%clearvars;

%% load data
% fileName = 'pancrea3_ch0_stack0000_488nm_0000000msec_0005772506msecAbs.tif';
% filePath = fullfile(userpath, 'ca_test_data', 'pancrea3_z1um_40ms', fileName);
fileName = 'proc_cell4a_zp24um_405.tif';
filePath = fullfile(userpath, fileName);
tic;
I = tiff.imread(filePath, true);
t = toc;
fprintf('data loaded in %.2fms\n', t*1e3);

%% show the data
h = gui.VolView( ...
    'Title', fileName, ...
    'Data', I, ...
    'VoxelSize', [0.103, 0.103, 0.24] ...
);

% set(gcf, 'WindowButtonDownFcn', @mouseDown);
% 
% function mouseDown(object, eventdata)
%     C = get(gca, 'CurrentPoint');
%     disp(['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
% end