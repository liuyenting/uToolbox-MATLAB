function ts = timestamp(format)
%TIMESTAMP Get timestamp.
%   
%   TS = TIMESTAMP create a string based on local time.
%   TS = TIMESTAMP(FORMAT) can further control the format of the output.
%   There are three types of plausible output:
%       - 'DateTime':   'yyyy-MM-dd-HH-mm-ss'
%       - 'Date':       'yyyy-MM-dd'
%       - 'Time':       'HH-mm-ss'
%
%   Note
%   ----
%   Time zone is currently locked with the system. This can be expanded by
%   passing additional parameter to the DATETIME function.
%
%   See also: DATETIME

if nargin == 0
    format = 'DateTime';
end

DATE_FORMAT = 'yyyy-MM-dd';
TIME_FORMAT = 'HH-mm-ss';
if strcmp(format, 'Time')
    format = TIME_FORMAT;
elseif strcmp(format, 'DateTime')
    format = [DATE_FORMAT, '-', TIME_FORMAT];
elseif strcmp(foramt, 'Date')
    format = DATE_FORMAT;
end

t = datetime('now', 'TimeZone', 'local');
ts = datestr(t, format);

end

