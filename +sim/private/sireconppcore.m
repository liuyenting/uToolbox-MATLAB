function J = sireconppcore(I, M, kp, parms)
%SIRECONPPCORE Summary of this function goes here
%   
%   TBA

%% parameters
imSz = size(I);
% remove the phases
imSz = imSz(1:2);

nOri = parms.Orientations;
nPhase = parms.Phases;
padSz = parms.PadSize;

pxSz = parms.PixelSize;

% interpolated size
rSz = parms.RetrievalInterpRatio*imSz;

%% pre-allocate
TF = zeros([imSz, nPhase], 'single');
k0 = zeros([nPhase, nOri], 'single');

TF(:, :, 1) = parms.TransFunc(:, :, 1); 

%% find initial phase
if parms.Debug
    dispFlag = 'iter-detailed';
else
    dispFlag = 'none';
end
for iOri = 1:nOri
    % convert to frequency space
    D = fftshift(fft2(ifftshift(I(:, :, :, iOri))));
    
    %% retrieve domains
    % flatten the array
    D = reshape(D, [prod(imSz), nPhase]);
    % solve the matrix
    D = (M \ D')';
    % reshape back to original image size
    D = reshape(D, [imSz, nPhase]);
    
    for iPhase = 2:nPhase
        O0 = parms.TransFunc(:, :, 1);
        D0 = D(:, :, 1);
        Om = parms.TransFunc(:, :, iPhase);
        Dm = D(:, :, iPhase);
        
        % grids for the relative phase shift matrix
        [vx, vy] = meshgrid(1:imSz(1), 1:imSz(2));
        midpt = floor(imSz/2)+1;
        vx = (vx - midpt(1)) / imSz(1);
        vy = (vy - midpt(2)) / imSz(2);
        
        %% shift the result in real-space
        kx = kp(1, iPhase-1, iOri);
        ky = kp(2, iPhase-1, iOri);
        
        shift = exp(-1i * 2*pi * (kx*vx + ky*vy));
        
        % shift the transfer function
        Om = fftshift(ifft2(ifftshift(Om)));
        Om = Om .* shift;
        Om = fftshift(fft2(ifftshift(Om)));
        
        % shift the domain
        Dm = fftshift(ifft2(ifftshift(Dm)));
        Dm = Dm .* shift;
        Dm = fftshift(fft2(ifftshift(Dm)));
        
        %% apply the compensation factors
        D0 = D0 .* Om;
        Dm = Dm .* O0;
        
        %% generat the mask
        % grid
        [vx, vy] = meshgrid(1:imSz(1), 1:imSz(2));
        midpt = floor(imSz/2)+1;
        vx = vx - midpt(1);
        vy = vy - midpt(2);
        % distance map
        DM0 = hypot(vx, vy);
        DMm = hypot(vx + kx, vy + ky);
        
        % cutoff frequency
        f = 2*parms.NA / parms.Wavelength;
        % cut-off radius in pixel
        r = f * (imSz.*pxSz);
        r = min(r(1));
        
        mask = zeros(imSz, 'single');
        mask((DM0 < r) & (DMm < r)) = 1;
       
        D0 = D0 .* mask;
        Dm = Dm .* mask;
        
        %% linear regression for the coefficient
        options = optimoptions(@lsqnonlin, ...
                               'Display', dispFlag, ...
                               'FiniteDifferenceType', 'central', ...
                               'FiniteDifferenceStepSize', [1e-2, 1e-3]);
        st = lsqnonlin(@(s)(errfunc(s, D0, Dm)), [1, 0], [], [], options);
        fprintf('\ts = %.4f * exp(1i * %.4f)\n', st(1), st(2));
        
        %% save the constant coefficient
        k0(iPhase, iOri) = st(1) * exp(1i * st(2));
    end
end

%% integrate the coefficients
% average the coefficients along the orientation
k0 = mean(k0, 2);

% integrate the initial phase shifts
TF(:, :, 2:end) = TF(:, :, 2:end) .* reshape(k0(2:end), 1, 1, []);

%% multiply the transfer function
% the original apodization function
A = genapomask(imSz, 0.8);

% buffer space for the Wiener filter functions
C = zeros(imSz, 'single');
% temporary padded result
Tp = zeros(rSz, 'single');

% merged result
R = zeros(rSz, 'single');

% filter function
for iOri = 1:nOri
    for iPhase = 1:nPhase
        %% variables
        if iPhase > 1
            kx = kp(1, iPhase-1, iOri);
            ky = kp(2, iPhase-1, iOri);
        else
            kx = 0; 
            ky = 0;
        end
        Om = TF(:, :, iPhase);
        
        % grids for the relative phase shift matrix
        [vx, vy] = meshgrid(1:imSz(1), 1:imSz(2));
        midpt = floor(imSz/2)+1;
        vx = (vx - midpt(1)) / imSz(1);
        vy = (vy - midpt(2)) / imSz(2);
        
        %% denominator
        for iiOri = 1:nOri
            for iiPhase = 1:nPhase
                fprintf('d=%d, m=%d, d''=%d, m''=%d\n', iOri, iPhase, iiOri, iiPhase);
                if iiPhase > 1
                    kix = kp(1, iiPhase-1, iiOri);
                    kiy = kp(2, iiPhase-1, iiOri);
                else
                    kix = 0; 
                    kiy = 0;
                end
                
                % calcualte the phase shift
                kdx = kix - kx;
                kdy = kiy - ky;
                
                shift = exp(-1i * 2*pi * (kdx*vx + kdy*vy));

                % shift the transfer function
                Oms = fftshift(ifft2(ifftshift(Om)));
                Oms = Oms .* shift;
                Oms = fftshift(fft2(ifftshift(Oms)));

                % squared absolute value
                Oms = abs(Oms).^2;
                % sum the result
                C = C + Oms;
            end
        end
        % apply the Wiener coefficient
        C = C + (parms.WienerConstant)^2;       
    end
end

for iOri = 1:nOri
    for iPhase = 1:nPhase
        %% variables
        if iPhase > 1
            kx = kp(1, iPhase-1, iOri);
            ky = kp(2, iPhase-1, iOri);
        else
            kx = 0; 
            ky = 0;
        end
        Om = TF(:, :, iPhase);
        Dm = D(:, :, iPhase);
        
        % grids for the relative phase shift matrix
        [vx, vy] = meshgrid(1:imSz(1), 1:imSz(2));
        midpt = floor(imSz/2)+1;
        vx = (vx - midpt(1)) / imSz(1);
        vy = (vy - midpt(2)) / imSz(2);
        
        %% nominator
        shift = exp(-1i * 2*pi * -(kx*vx + ky*vy));

        % shift the transfer function
        As = fftshift(ifft2(ifftshift(A)));
        As = As .* shift;
        As = fftshift(fft2(ifftshift(As)));

        % generate the nominator (negative sign)
        N = conj(Om) .* As;
        
        % apply the transfer function
        T = (N./C) .* Dm;
        
        % pad the result to center
        li = floor((rSz-imSz)/2)+1;
        ui = li+imSz-1;
        Tp(li(1):ui(1), li(2):ui(2)) = T;
        
        % grids for the relative phase shift matrix
        [vx, vy] = meshgrid(1:rSz(1), 1:rSz(2));
        midpt = floor(rSz/2)+1;
        vx = (vx - midpt(1)) / rSz(1);
        vy = (vy - midpt(2)) / rSz(2);
        
        shift = exp(-1i * 2*pi * -(kx*vx + ky*vy));
        
        % shift in real space with complex gradient
        Tp = fftshift(ifft2(ifftshift(Tp)));
        Tp = Tp .* shift;
        Tp = fftshift(fft2(ifftshift(Tp)));
        
        %% add the result
        R = R + Tp;
    end
end

%% revert back to real space and show the result
R = fftshift(ifft2(ifftshift(R)));

figure('Name', 'Result', 'NumberTitle', 'off');
imagesc(abs(R));
    axis image;
    colormap(gray);
disp('DONE');




%% pre-allocate
% buffer space for results from the frequency domain, each for the original
% image and the padded image
F = zeros([imSz, nPhase], 'single');
Fp = zeros([imSz, nPhase, nOri], 'single');

% initial phase 
popt = zeros([(nPhase-1)/2, nOri], 'single');

% buffer space for the relative matrix (single orientation only)
pr = zeros([rSz, nPhase-1, nOri], 'single');
% grids for the relative phase shift matrix
[vx, vy] = meshgrid(1:rSz(1), 1:rSz(2));
midpt = floor(rSz/2)+1;
vx = vx - midpt(1);
vy = vy - midpt(2);

%% process
for iOri = 1:nOri
    fprintf('.. o = %d\n', iOri);
    
    %% create relative phase shift matrix
    % revert from position to shift
    midpt = floor(imSz/2)+1;
    shift = kp(:, :, iOri) - midpt.';
    % the ratio in current upsampled dimension
    shift = bsxfun(@rdivide, shift, rSz.');
    % calculate shift in unit spatial frequency, the result is negated
    % since we are trying to shift it back to where it should be
    shift = (2*pi) * (-shift);
    
    for iPhase = 1:nPhase-1
        % fill-in the distance matrix with phase shifts (in unit spatial
        % frequency)
        D = vx*shift(1, iPhase) + vy*shift(2, iPhase);
        % convert to imaginary part in order to apply shift in R-space
        pr(:, :, iPhase, iOri) = exp(1i * D);
    end
    
    for iPhase = 1:nPhase
        % extract volume
        T = I(:, :, iPhase, iOri);
        
        % pad the surrounding sides
        T = padarray(T, [padSz, padSz], 0, 'both');
        % RL deconvolution
        T = deconvlucy(T, parms.PSF, parms.PreDeconv);
        % crop the result
        T = T(padSz+1:end-padSz, padSz+1:end-padSz);
        
        % perform Fourier transform without the padded region
        F(:, :, iPhase) = fftshift(fft2(ifftshift(T)));
    end
    
    %% retrieve domains
    % flatten the array
    F = reshape(F, [prod(imSz), iPhase]);
    % solve the matrix
    F = (M \ F')';
    % reshape back to original image size
    F = reshape(F, [imSz, iPhase]);
    
    %% find phases
%     % multiply apodization function
%     F = F .* A;
%     % upsampling to perform FT interpolation in real space
%     li = floor((rSz-imSz)/2)+1;
%     ui = li+imSz-1;
%     Fp(li(1):ui(1), li(2):ui(2), :, iOri) = F;
    Fp(:, :, :, iOri) = F;
    
%     %% test run
%     nt = 20;
%     pt = linspace(0, 2*pi, nt);
%     ct = zeros([nt, 1]);
%     for t = 1:nt
%         p = pt(t); 
%         ct(t) = costfunc(Fref, Fp(:, :, 4:5), p, pr(:, :, 3:4));
%         fprintf('t = %d, p1 = %f, c = %f\n', t, p, ct(t));
%         pause(2);
%     end
%     figure('Name', 'Cost Function t-Plot', 'NumberTitle', 'off');
%     plot(ct);
%         xlabel('Phase Shift');
    
    %% search the optimal inital phase
    % unit spatial frequency
    lim = (2*pi) ./ imSz;
    
    options = optimoptions( ...
        'fmincon', ...
        'FiniteDifferenceStepSize', max(lim), ...
        'StepTolerance', 1e-2, ...
        'Display', 'notify-detailed' ...
    );

    p0 = fmincon( ...
        @(x) costfunc(F(:, :, 1), F(:, :, 2:end), imSz, rSz, x, pr), ...
        [0, 0], ...
        [], [], [], [], ...
        [-pi, -pi], [pi, pi], ...
        [], ...
        options ...
    );

    % round p0 to nearest multiple of unit spatial frequency  
    p0 = round(p0./lim) .* lim;
    fprintf('.... m1=%f, m2=%f\n', p0(1), p0(2));
    
    % save the optimal initial phase shift
    popt(:, iOri) = p0;
end

%% the actual reconstruction
%TODO use generalized Weiner filter
% % sum all the orientations and phases
% Fopt = reshape(Fopt, [rSz, nOri*nPhase]);
% J = sum(Fopt, 3);
% J = fftshift(ifft2(ifftshift(J), 'symmetric'));

%TODO temporary smaller relative phase shift matrix

% buffer space for the relative matrix (single orientation only)
prs = zeros([imSz, nPhase-1, nOri], 'single');
% grids for the relative phase shift matrix
[vx, vy] = meshgrid(1:imSz(1), 1:imSz(2));
midpt = floor(imSz/2)+1;
vx = vx - midpt(1);
vy = vy - midpt(2);

for iOri = 1:nOri
    % revert from position to shift
    midpt = floor(imSz/2)+1;
    shift = kp(:, :, iOri) - midpt.';
    % the ratio in current upsampled dimension
    shift = bsxfun(@rdivide, shift, imSz.');
    % calculate shift in unit spatial frequency, the result is negated
    % since we are trying to shift it back to where it should be
    shift = (2*pi) * (-shift);
    
    for iPhase = 1:nPhase-1
        % fill-in the distance matrix with phase shifts (in unit spatial
        % frequency)
        D = vx*shift(1, iPhase) + vy*shift(2, iPhase);
        % convert to imaginary part in order to apply shift in R-space
        prs(:, :, iPhase, iOri) = exp(1i * D);
    end
end

J = wnrrecon(F, imSz, popt, prs, pr, parms);

end

function err = costfunc(Fref, F, imSz, rSz, p0, pr)
%COSTFUNC Cost function to minimize for the phase retrieval algorithm.
%
%   Fref: reference frequency plane, m_0
%     Fp: frequency plane, m_i
%     p0: initial phase shift
%     pr: relative phsae shift, determined by kp

% persistent h;
% 
% if isempty(h) || ~isvalid(h)
%     h = figure('Name', 'Phase Retrieval', 'NumberTitle', 'off');
% end

profile resume;

% interleave the phases since we now have m_i^- and m_i^+
p0 = exp(1i * [+p0; -p0]);
% flatten the array for linear duplication
p0 = p0(:);
np = length(p0);

F = F .* reshape(p0, [1, 1, np]);

Fp = zeros([rSz, np], 'single');
% upsampling to perform FT interpolation in real space
li = floor((rSz-imSz)/2)+1;
ui = li+imSz-1;
Fp(li(1):ui(1), li(2):ui(2), :) = F;

% back to time domain
Fp = fftshift(ifft2(ifftshift(Fp)));

% add relative phase shift deduced from kp values (imaginary number in the
% time domain)
Fp = Fp .* pr;

Fp = fftshift(fft2(ifftshift(Fp)));

% % sum the result to evaluate performance
% S = Fref + sum(Fp, 3);
% 
% R = fftshift(ifft2(ifftshift(S), 'symmetric'));
% figure(h);
% imagesc(R);
%     axis image;
%     colormap(gray);
% drawnow;

Frefp = zeros(rSz, 'single');
% upsampling to perform FT interpolation in real space
li = floor((rSz-imSz)/2)+1;
ui = li+imSz-1;
Frefp(li(1):ui(1), li(2):ui(2)) = Fref;

% error
err = abs(Frefp - Fp);
err = sum(err(:)) / np;

% output is required to be double instead of single
err = double(err);

profile off;

end

function err = errfunc(s, a, b)
%ERRFUNC Cost function for the non-linear fitting of phi.
%
%   TBA

err = b.*(s(1) * exp(1i*s(2))) - a;
% return type is required to be double
err = double(err);

end
