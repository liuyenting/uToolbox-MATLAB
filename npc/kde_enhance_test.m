close all;
clearvars;

filePath = 'cell4a_zp24um_405.tif';

[I, sz] = tiff.imread(filePath);
J = zeros(sz, 'uint16');

sz

for iz = 1:sz(3)
    disp(iz);
    Ip = I(:, :, iz);
    [Jp, ~] = wiener2(Ip, [5, 5]);
    J(:, :, iz) = uint16(Jp);
end

tiff.imsave(J, ['proc_', filePath], true);