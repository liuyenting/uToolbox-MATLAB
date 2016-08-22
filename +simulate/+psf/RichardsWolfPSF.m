classdef RichardsWolfPSF < AbstractPSF
    %RICHARDSWOLFPSF Generate a PSF using the Richards & Wolf model.
    
    properties (Access = private)
        rGrid;
    end
    
    methods (Access = public)
        function obj = RichardsWolfPSF
            obj = obj@AbstractPSF('Richards & Wolf');
        end
        
        function obj = generate(obj)
            nz = obj.nz;
            
            % Find out the center.
            midz = (nz-1)/2;
            
            % Generate the calculation grid along radial direction.
            obj.generatergrid();

            % Generate the PSF plane by plane.
            for iz = 1:obj.nz
                % Convert to real world unit.
                rz = obj.axialRes * (iz - midz);
                
                data = obj.generatelayer(rz);
                obj.data.setplane(iz, data);
                disp(iz);
            end
        end
    end
    
    methods (Access = private)
        function obj = generatergrid(obj)
            % Generate the radial grid for calculation.
            
            % Find out the center.
            sz = [obj.nx, obj.ny];
            cp = (sz-1)/2;
            
            % Calculate the maximum radius in pixel.
            maxRadius = floor(sqrt( sum((sz-cp).^2) )) + 1;
            
            % Expand the grid size according to the over sampling ratio.
            obj.rGrid = 1:(maxRadius*obj.OVER_SMPL_RATIO);
            obj.rGrid = (obj.rGrid-1) / obj.OVER_SMPL_RATIO;
            
            % Convert to real world unit in micron.
            obj.rGrid = obj.rGrid * obj.radialRes;
        end
        
        function data = generatelayer(obj, z)
            % Generate PSF of plane with defocus distance Z in micron.

            % Calculate the Kirchhoff's diffraction value along the grid.
            kdIntegral = KirchhoffDiffraction(obj.NA, obj.lambda, z, obj.ni);
            hLut = arrayfun(@kdIntegral.calculate, obj.rGrid);
            
            % Convert to polar coordinates.
            sz = [obj.nx, obj.ny];
            cp = (sz-1)/2;
            [xMat, yMat] = meshgrid(1:sz(1), 1:sz(2));
            [~, rMat] = cart2pol(xMat-cp(1), yMat-cp(2));
            
            % Convert to real world unit.
            rMat = rMat.' * obj.radialRes;
            
            % Fill the data.
            data = interp1(obj.rGrid, hLut, rMat, 'linear', 'extrap');
        end
    end
    
end

