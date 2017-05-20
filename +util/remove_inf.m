%% open the file
inPath = fullfile(userpath, 'FRC_05202017', 'Exposuretime_dependent(sigle layer)', '75msec_result_driftcorrected.csv');
fprintf('path = "%s"\n', inPath);

outPath = util.chfext(inPath, 'red.csv');

fin = fopen(inPath, 'r');
fout = fopen(outPath, 'w+');

%% process line by line
tline = fgets(fin);
skipped = 0;
while ischar(tline)
    if isempty(regexp(tline,'[¡Û]', 'once'))
        fwrite(fout, tline);
    else
        skipped = skipped+1;
    end
    tline = fgets(fin);
end

fprintf('skipped %d lines\n', skipped);

%% cleanup
fclose(fin);
fclose(fout);