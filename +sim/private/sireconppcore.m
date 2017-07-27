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

hMask = figure('Name', 'Mask', 'NumberTitle', 'off');

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
        
        %% shift the result in real-space
        kx = kp(1, iPhase-1, iOri);
        kx = round(kx); 
        ky = kp(2, iPhase-1, iOri);
        ky = round(ky);
        
        % translated frame
        li = 1 + [kx, ky];
        li(li < 1) = 1;
        ui = imSz + [kx, ky];
        ui(ui > imSz) = imSz(1);
        % original frame
        l0i = li - [kx, ky];
        u0i = ui - [kx, ky];
        
        Om = zeros(imSz, 'single');
        Om(li(2):ui(2), li(1):ui(1)) = parms.TransFunc(l0i(2):u0i(2), l0i(1):u0i(1), iPhase);
        
        Dm = zeros(imSz, 'single');
        Dm(li(2):ui(2), li(1):ui(1)) = D(l0i(2):u0i(2), l0i(1):u0i(1), iPhase);

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
        DMm = hypot(vx - kx, vy - ky);
        
        % cutoff frequency
        f = 2*parms.NA / parms.Wavelength;
        % cut-off radius in pixel
        r = f * (imSz.*pxSz);
        r = min(r(1));
        
        mask = zeros(imSz, 'single');
        mask((DM0 < r) & (DMm < r)) = 1;
        
        % preview mask
        figure(hMask);
        % generate title string
        m = floor(iPhase/2);
        if iPhase > 1
            if mod(iPhase, 2) == 0
                s = '^-';
            else
                s = '^+';
            end
        else
            s = '';
        end
        t = sprintf('d_%d, m_%d%s', iOri, m, s);
        subplot(nOri, nPhase-1, (iOri-1)*(nPhase-1)+iPhase-1);
        imagesc(mask);
            axis image;
            title(t);
        drawnow;
       
        D0 = D0 .* mask;
        Dm = Dm .* mask;
        
        %% linear regression for the coefficient
        % nominator
        nom = conj(Dm) .* D0;
        nom = sum(nom(:));
        % denominator
        den = abs(Dm).^2;
        den = sum(den(:));
        % constant complex coefficient
        s = nom / den;
        fprintf('\ts = %.4f * exp(1i * %.4f)\n', real(s), imag(s));
        
        %% save the constant
        k0(iPhase, iOri) = s;
    end
end

%% integrate the coefficients
% average the coefficients along the orientation
k0 = mean(k0, 2);

TF(:, :, 1) = parms.TransFunc(:, :, 1); 
% integrate the initial phase shifts
TF(:, :, 2:end) = parms.TransFunc(:, :, 2:end) .* reshape(k0(2:end), 1, 1, []);

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
        
        %% denominator
        for iiOri = 1:nOri
            for iiPhase = 1:nPhase
                if iiPhase > 1
                    kix = kp(1, iiPhase-1, iiOri);
                    kiy = kp(2, iiPhase-1, iiOri);
                else
                    kix = 0; 
                    kiy = 0;
                end
                
                % calcualte the phase shift
                kix = kix - kx;
                kix = round(kix);
                kiy = kiy - ky;
                kiy = round(kiy);

                li = 1 + [kix, kiy];
                li(li < 1) = 1;
                ui = imSz + [kix, kiy];
                ui(ui > imSz) = imSz(1);

                l0i = li - [kix, kiy];
                u0i = ui - [kix, kiy];
                
                Oms = zeros(imSz, 'single');
                Oms(li(2):ui(2), li(1):ui(1)) = Om(l0i(2):u0i(2), l0i(1):u0i(1));

                % squared absolute value
                Oms = abs(Oms).^2;
                % sum the result
                C = C + Oms;
            end
        end
        % apply the Wiener coefficient
        C = C + (parms.WienerConstant)^2;       

        Om = TF(:, :, iPhase);
        Dm = D(:, :, iPhase);
        
        kx = round(kx); 
        ky = round(ky);

        li = 1 + [kx, ky];
        li(li < 1) = 1;
        ui = imSz + [kx, ky];
        ui(ui > imSz) = imSz(1);

        l0i = li - [kx, ky];
        u0i = ui - [kx, ky];
        
        %% nominator
        As = zeros(imSz, 'single');
        As(li(2):ui(2), li(1):ui(1)) = A(l0i(2):u0i(2), l0i(1):u0i(1));

        % generate the nominator (negative sign)
        N = conj(Om) .* As;
        
        % apply the transfer function
%         T = (N./C) .* Dm;
        T = Dm;
        if iPhase > 1
            t = k0(iPhase, iOri);
            T = T .* exp(1i * atan(imag(t)/real(t)));
        end
        
        % pad the result to center
        li = floor((rSz-imSz)/2)+1;
        ui = li+imSz-1;
        Tp = zeros(rSz, 'single');
        Tp(li(1):ui(1), li(2):ui(2)) = T;

        % shift in real space with complex gradient
        Tp = fftshift(ifft2(ifftshift(Tp)));
        
        % grids for the relative phase shift matrix
        [vx, vy] = meshgrid(1:rSz(1), 1:rSz(2));
        midpt = floor(rSz/2)+1;
        vx = (vx - midpt(1)) / rSz(1);
        vy = (vy - midpt(2)) / rSz(2);
        
        if iPhase > 1
            kx = kp(1, iPhase-1, iOri);
            ky = kp(2, iPhase-1, iOri);
        else
            kx = 0; 
            ky = 0;
        end
        shift = exp(-1i * 2*pi * (kx*vx + ky*vy));
        Tp = Tp .* shift;
        
        Tp = fftshift(fft2(ifftshift(Tp)));
        
        %% add the result
        R = R + Tp;
    end
end

%% revert back to real space and show the result
figure('Name', 'Result', 'NumberTitle', 'off');
subplot(1, 2, 1);
imagesc(abs(R).^0.1);
    axis image;
    title('Frequency Domain');
    
J = fftshift(ifft2(ifftshift(R)));

subplot(1, 2, 2);
imagesc(abs(J));
    axis image;
    title('Time Domain');
    
disp('DONE');

end
