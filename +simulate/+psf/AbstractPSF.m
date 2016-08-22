classdef (Abstract) AbstractPSF < handle
    %ABSTRACTPSF Abstract class for the point spread functions.
    
    properties (GetAccess = public, SetAccess = protected)
        data;       % Generated PSF data.
    end
    
    properties (Access = protected)
        type;       % Type of the PSF model.
        
        NA;         % Numerical aperture.
        lambda;     % Wavelength in nanometer.
        ni;         % Refraction index.
        
        axialRes;   % Resolution along axial direction in micron.
        radialRes;  % Resolution along radial direction in micron.
    end
    
    properties (Access = protected)
        nx;         % X dimension of the PSF in pixel.
        ny;         % Y dimension of the PSF in pixel.
        nz;         % Z dimension of the PSF in pixel.
    end
    
    properties (Access = protected, Constant = true)
        OVER_SMPL_RATIO = 10;     % Over sampling ratio.
        PLOT_STEPS = 256;         % Steps of the generated plots.
        AXIALTICK_DENSITY = 2;    % AXIALTICK_DENSITY lambdas per X tick.
        RADIALTICK_DENSITY = 1/2; % RADIALTICK_DENSITY labmdas per Y tick.
    end
    
    methods (Access = public)
        function obj = AbstractPSF(type)
            obj.type = type;
        end
        
        function obj = setopticsparam(obj, NA, lambda, ni)
            % Set the optical parameters.
            
            obj.NA = NA;
            obj.lambda = lambda;
            obj.ni = ni;
        end
        
        function obj = setresparam(obj, axialRes, radialRes)
            % Set the resolution parameters.
            
            obj.axialRes = axialRes;
            obj.radialRes = radialRes;
        end
        
        function obj = setoutputparam(obj, dim)
            % Set the output dimension. Using this function will wipe
            % previous calculation result!
            
            % Save the dimension size to local storage.
            obj.nx = dim(1);
            obj.ny = dim(2);
            obj.nz = dim(3);
            
            % Create the storage object.
            if isempty(obj.data)
                warning('PSFGenerator:AbstractPSF:setoutputparam', ...
                        'Volume object overwrited.');
            end
            obj.data = VolumeObj(dim);
        end
    end
    
    methods (Abstract)
        obj = generate(obj);
    end
    
    methods
        function plotortho(obj, dim)
            dim = upper(dim);
            
            slice = obj.data.getdata();
            
            % Parse the parameter. 
            switch dim
                case 'XY'
                    cp = floor((obj.nz-1)/2) + 1;
                    slice = slice(:, :, cp);
                    
                    xtickdensity = obj.RADIALTICK_DENSITY;
                    ytickdensity = obj.RADIALTICK_DENSITY;
                    
                    xres = obj.radialRes;
                    yres = obj.radialRes;
                    
                    xmax = obj.ny;
                    ymax = obj.nx;
                case 'XZ'
                    cp = floor((obj.ny-1)/2) + 1;
                    slice = slice(:, cp, :);
                    
                    xtickdensity = obj.AXIALTICK_DENSITY;
                    ytickdensity = obj.RADIALTICK_DENSITY;
                    
                    xres = obj.axialRes;
                    yres = obj.radialRes;
                    
                    xmax = obj.nz;
                    ymax = obj.nx;
                case 'YZ'
                    cp = floor((obj.nx-1)/2) + 1;
                    slice = slice(cp, :, :);
                    
                    xtickdensity = obj.AXIALTICK_DENSITY;
                    ytickdensity = obj.RADIALTICK_DENSITY;
                    
                    xres = obj.axialRes;
                    yres = obj.radialRes;
                    
                    xmax = obj.nz;
                    ymax = obj.ny;
                otherwise
                    error('PSFGenerator:AbstractPSf:plotortho', ...
                          'Unknown orthoslice dimension.');
            end
            xcp = floor((xmax-1)/2) + 1;
            ycp = floor((ymax-1)/2) + 1;
            
            % Squeeze the singleton dimension.
            slice = squeeze(slice);
            % Normalize the slice to 256 steps.
            maxVal = max(slice(:));
            slice = (obj.PLOT_STEPS-1) * (slice/maxVal);
            
            % Plot in a new figure.
            figureTitle = [obj.type, ' Orthoslice'];
            figure('Name', figureTitle, 'NumberTitle', 'off');
            image(slice);
            colormap(jet(obj.PLOT_STEPS));
            % Set the aspect ratio.
            pbaspect([xres*xmax, yres*ymax, 1]);
            
            % Recalculate ticks position.
            xtickInt = xtickdensity * (obj.lambda*1e-3) / xres;
            xtickLo = xcp:-xtickInt:1;
            xtickHi = xcp:xtickInt:xmax;
            xtick = [fliplr(xtickLo), xtickHi(2:end)];
            
            ytickInt = ytickdensity * (obj.lambda*1e-3) / yres;
            ytickLo = ycp:-ytickInt:1;
            ytickHi = ycp:ytickInt:ymax;
            ytick = [fliplr(ytickLo), ytickHi(2:end)];
            
            % Convert the tick labels from pixel to micron.
            xticklabel = (xtick-xcp) * xres;
            yticklabel = (ytick-ycp) * yres;
            
            % Set the axes.
            set(gca, 'XTick', xtick);
            set(gca, 'XTickLabel', sprintf('%.3f\n', xticklabel));
            
            set(gca, 'YTick', ytick);
            set(gca, 'YTickLabel', sprintf('%.3f\n', yticklabel));
  
            % Overlay figure object with different tick color.
            a = gca;
            b = copyobj(a, gcf);
            set(b, 'Xcolor', 'white', 'YColor', 'white', ...
                   'XTickLabel', [], 'YTickLabel', []);
               
            % Set the view title.
            % Note: Delayed till here to avoid duplication.
            title(dim);
            xlabel([dim(2), ' (\mum)']);
            ylabel([dim(1), ' (\mum)']);
        end
    end
    
end

