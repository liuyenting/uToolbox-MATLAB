function dispkp(kp, parms)
%DISPKP Display pattern wave vector results in a table.
%   
%   DISPKP(KP, PARMS) dumps the wave vectors in a formulated table --
%   orientations as columns and each term as a row, with coordinates
%   written in pairs.

%% parameters
nOri = parms.Orientations;
nPhase = parms.Phases;

%% generate table
% generate column and row labels
colname = cell([1, nOri]);
for iOri = 1:nOri
    colname{iOri} = sprintf('Orientation%d', iOri);
end
rowname = cell([1, nPhase-1]);
for iPhase = 2:2:nPhase
    i = iPhase/2;
    rowname{iPhase-1} = sprintf('m%d-', i);
    rowname{iPhase} = sprintf('m%d+', i);
end

% generate (x, y) coordinate pair
kpstr = cell([nPhase-1, nOri]);
for iOri = 1:nOri
    for iPhase = 1:nPhase-1
        kpstr{iPhase, iOri} = sprintf( ...
            '(%.4f, %.4f)', kp(1, iPhase, iOri), kp(2, iPhase, iOri) ...
        );
    end
end

% create the table and print-out
result = array2table(kpstr, ...
                     'VariableNames', colname, 'RowNames', rowname);

%% print the result
fprintf('\n');
disp(result);

end
