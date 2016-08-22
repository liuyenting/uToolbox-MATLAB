classdef KirchhoffDiffraction
    %KIRCHHOFFDIFFRACTION Calculate the Kirchhoff's diffraction formula.
    
    properties (Access = private)
        NA;     % Numerical aperture.
        lambda; % Wavelength in nanometer.
        z;      % Distance from the focus on the optical axis in micron.
        ni;     % Refraction index (unused).
    end
    
    methods (Access = public)
        function obj = KirchhoffDiffraction(NA, lambda, z, ni)
            obj.NA = NA;
            obj.lambda = lambda;
            obj.z = z;
            obj.ni = ni;
        end
        
        function I = calculate(obj, r)
            % Calculate the intensity of distance R (micron) from the 
            % optical axis along the raidal direction.

            E = integral(@(rho)obj.integrand(rho, r), 0, 1);
            I = abs(E)^2;
        end
    end
    
    methods (Access = private)
        function i = integrand(obj, rho, r)
            % The integrand of the Kirchhoff's diffraction formula.
            
            % Wave number.
            k0 = (2*pi) / obj.lambda;
            %k = obj.ni * k0;
            k = k0;
            
            % Bessel function.
            x = k * (obj.NA/obj.ni) * (r*1e3) * rho;
            J0 = besselj(0, x);
            
            % Optical path difference.
            OPD = -1/2 * (obj.NA/obj.ni)^2 * (obj.z*1e3) * rho.^2;
            
            % Complex polar coordinate.
            theta = k * OPD;
            i = J0 .* exp(1i * theta) .* rho;
        end
    end
    
end

