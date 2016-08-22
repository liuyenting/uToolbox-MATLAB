classdef Shear < util.cppbridge.CppBridge
    %SHEAR Shearing wrapper class.

    properties (Constant, Access = private, Hidden = true)
        MSG_ID = 'postproc:shear:matlab';
        MEX_NAME = 'shearmex';
    end

    properties (Access = private)
        % Original intensity scale of the image.
        intensity;
    end

    methods
        function this = Shear()
            this@util.cppbridge.CppBridge(mfilename('fullpath'), ...
                                          postproc.Shear.MEX_NAME);
        end

        function this = setacqparam(this, param)
            %SETPARAM Set the acquisition paramter.

            % All the fields should be of type SINGLE.
            param = structfun(@single, param, 'UniformOutput', false);

            this.mexHandle('setacqparam', this.objectHandle, ...
                           param);
        end

        function this = loadstack(this, im)
            %IMLOAD Load the image stack to be processed.

            % Image should be a 3-D array of type UINT16.
            if ndims(im) ~= 3
                error(postproc.Shear.MSG_ID, ...
                      'Image should be a stack.');
            elseif ~isa(im, 'uint16')
                error(postproc.Shear.MSG_ID, ...
                      'Image should be of type UINT16.');
            end

            %% Upload to device.
            im = gpuArray(im);

            %% Probe the original intensity range.
            this.intensity.min = single(min(im(:)));
            this.intensity.max = single(max(im(:)));

            %% Permute from X-Y-Z to Y-Z-X.
            im = permute(im, [2, 3, 1]);

            this.mexHandle('loadstack', this.objectHandle, im);
        end

        function this = execute(this)
            %EXECUTE Invoke the shear kernel for the loaded image.

            this.mexHandle('execute', this.objectHandle);
        end

        function im = retrieveresult(this)
            %GETRESULT Retrieve result from the post processing pipeline.

            im = this.mexHandle('retrieveresult', this.objectHandle);

            %% Adjust the intensity.
            im = this.scaleIntensity(im);

            %% Permute from Y-Z-X to X-Y-Z.
            im = permute(im, [3, 1, 2]);

            %% Pull data from device.
            im = gather(im);
        end
    end

    methods (Access = private)
        function J = scaleIntensity(this, I)
            %SCALEINTENSITY Scale the intensity of I to the range of
            %THIS.INTENSITY.

            mMin = min(I(:));
            mMax = max(I(:));

            oMin = this.intensity.min;
            oMax = this.intensity.max;

            ratio = (oMax-oMin) / (mMax-mMin);
            J = ratio*(I - mMin) + oMin;
        end
    end

    methods (Static, Access = protected)
        function hdl = getfunchandle()
            hdl = str2func(postproc.Shear.MEX_NAME);
        end
end
