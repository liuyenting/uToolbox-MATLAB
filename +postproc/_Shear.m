classdef Shear < handle
    %SHEAR Shearing wrapper class.
    
    properties (Constant, Access = private, Hidden = true)
        MATLAB_MSGID = 'postproc:shear:matlab';
        MEX_NAME = 'shearmex';
    end
    
    properties (SetAccess = immutable, GetAccess = private, Hidden = true)
        % Handle for the MEX function.
        mexHandle;
    end
    
    properties (Access = private, Hidden = true) 
        % Handle to the unerlying C++ class instance.
        objectHandle;
    end
    
    properties (Access = private)
        % Original intensity scale of the image.
        intensity;
    end

    methods
        function this = Shear()
            %CONSTRUCTOR Create a new C++ class instance.
            %
            %   Note: If the MEX file does not exist, we will have to
            %   compile it first.
            
            %% Test the MEX file availability.
            if ~postproc.Shear.ismexcompiled()
                postproc.Shear.compilemex();
            end
            
            %% Retrieve the function handle.
            this.mexHandle = str2func(postproc.Shear.MEX_NAME);
            
            %% Create the class instance.
            this.objectHandle = this.mexHandle('new');
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
                error(postproc.Shear.MATLAB_MSGID, ...
                      'Image should be a stack.');
            elseif ~isa(im, 'uint16') 
                error(postproc.Shear.MATLAB_MSGID, ...
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
        
        function delete(this)
            %DESTRUCTOR Call the desctructor of the created C++ class
            %instance.
            
            this.mexHandle('delete', this.objectHandle);
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
    
    methods (Static, Access = private)
        function b = ismexcompiled()
            %ISMEXCOMPILED Validate whether all the MEX files are compiled.

            % Search for the MEX file in path.
            b = exist(postproc.Shear.MEX_NAME, 'file') == 3;
        end
        
        function compilemex()
            %COMPILEMEX Compile all the necessary MEX files that will be
            %used by the wrapper class.

            %% Change directory into private folder.
            [currentDir, ~, ~] = fileparts(mfilename('fullpath'));
            targetDir = fullfile(currentDir, 'private');
            localDir = cd(targetDir);

            %% Switch back the directory when COMPILEMEX is exited.
            cdLocalObj = onCleanup(@() cd(localDir)); 

            %% Setup compile paramter.
            % CUDA path is hard coded for SDK v7.5, CUDA_INC_PATH is not
            % created for later CUDA SDKs, so we have to generate it by
            % ourselves.
            cudaIncPath = fullfile(getenv('CUDA_PATH_V7_5'), 'include');
            if exist(cudaIncPath, 'dir') ~= 7
                error(postproc.Shear.MATLAB_MSGID, ...
                      'CUDA_PATH is not properly configured.');
            else
                % Generate the include directive. 
                % No space shall place after the '-I' switch.
                cudaIncPath = ['-I''', cudaIncPath, ''''];
            end
            
            %% Compile the file.
            mexcuda('-v', ...                       % Verbose.
                    cudaIncPath, ...                % CUDA header path.
                    '-output', postproc.Shear.MEX_NAME, ...  % Filename.
                    '*.cpp', '*.cu');               % Source files.
        end
    end

end
