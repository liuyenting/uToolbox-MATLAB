function tokens = parspimfname(fname)
%PARSPIMFNAME Parse standard SPIM generated filename.
%   
%   TOKENS = PARSPIMFNAME(FNAME) parse the filename FNAME into a named
%   structure TOKENS, the following fields exist
%       - Project: Project name, used as prefix in most cases.
%       - Channel: Current acquisition channel, used with 'Wavelength.'
%       - StackNumber: Number of stacks in this time lapse.
%       - Wavelength: Wavelength in nanometer.
%       - RelativeTime: Relative timestamp in msec.
%       - AbsoluteTime: Absolute timestamp in msec.

expr = [
    '(?<Project>\w+)_' ...
    'ch(?<Channel>\d+)_' ...
    'stack(?<StackNumber>\d{4})_' ...
    '(?<Wavelength>\d+)nm_' ...
    '(?<RelativeTime>\d+)msec_' ...
    '(?<AbsoluteTime>\d+)msecAbs'
];
tokens = regexp(fname, expr, 'names');

end

