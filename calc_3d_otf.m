psf = 'data/08302017_SIM_test/RAWsingle bead_445_zp1um_step_p15um_50ms_ch0_stack0000_445nm_0000000msec_0000318072msecAbs.tif';
psf = fullfile(userpath, psf);

nori = 1;
nphases = 5;

Ipsf = io.fread(psf);
Ipsf = single(Ipsf);

Ipsf = sim.opmajor(Ipsf, nori, nphases);
volsz = size(Ipsf);
volsz = volsz(1:3);

Fpsf = zeros(size(Ipsf));
for i = 1:nphases
    Fpsf(:, :, :, i) = fftn(Ipsf(:, :, :, i));
end