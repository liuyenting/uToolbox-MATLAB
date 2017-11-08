function result = batchfrc
%BATCHFRC A simple wrapper utility to calculate FRC in an entire folder.
%
%   See also: RESOLUTION.RUNFRC

indir = uigetdir('C:\', 'Where are the coordinate files?');
% check whether the input directory exists
if exist(indir, 'dir') ~= 7
    error(generatemsgid('InvalidInDir'), ...
          'Input directory does not exists.');
end

% start the process
fds = fileDatastore(indir, ...
                    'ReadFcn', @resolution.runfrc, ...
                    'FileExtensions', '.csv');
result = readall(fds);

close all;

end