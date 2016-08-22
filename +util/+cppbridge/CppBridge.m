classdef CppBridge < handle
    %CPPBRIDGE Wrapper class for MATLAB/C++ object bridge.

    properties (Constant, Access = private, Hidden = true)
        MSGID = 'util:cppbridge:matlab';
        FORCE_COMPILE = false;
    end

    properties (SetAccess = protected, GetAccess = protected, Hidden = true)
        mexHandle;      % Handle for the MEX function.
        objectHandle;   % Handle to the unerlying C++ class instance.
        brdIncPath;     % Include directory for the bridge.
    end

    methods
        function this = CppBridge(clsPath, mexName)
            %CONSTRUCTOR Create a new C++ class instance.
            %   CPPBRIDGE(CLSPATH, MEXNAME) creates the bridge according
            %   the assigned class path CLSPATH and target MEX filename
            %   MEXNAME.
            %
            %   Note: If the MEX file does not exist, we will have to
            %   compile it first.
            
            %% Test the MEX file availability.
            [clsDir, ~, ~] = fileparts(clsPath);
            if ~util.cppbridge.CppBridge.ismexcompiled(clsDir, mexName)
                util.cppbridge.CppBridge.compilemex(clsDir, mexName);
            end

            %% Retrieve the function handle.
            this.mexHandle = this.getfunchandle();
            
            %% Create the class instance.
            this.objectHandle = this.mexHandle('new');
        end

        function delete(this)
            %DESTRUCTOR Call the desctructor of the created C++ class
            %instance.

            this.mexHandle('delete', this.objectHandle);
        end
    end

    methods (Static, Access = private)
        function b = ismexcompiled(clsDir, mexName)
            %ISMEXCOMPILED Validate whether all the MEX files are compiled.
            
            % Search for the MEX file in path.
            mexName = [mexName, '.', mexext];
            mexPath = fullfile(clsDir, 'private', mexName);
            b = exist(mexPath, 'file') && ...
                    ~util.cppbridge.CppBridge.FORCE_COMPILE;
        end

        function compilemex(clsDir, mexName)
            %COMPILEMEX Compile all the necessary MEX files that will be
            %used by the wrapper class.

            %% Extract target directory.
            targetDir = fullfile(clsDir, 'private');
            localDir = cd(targetDir);

            %% Switch back the directory when COMPILEMEX is exited.
            cdLocalObj = onCleanup(@() mexCleanup(localDir));
            function mexCleanup(d)
                % Delete all the object files.
                delete('*.obj');
                % Return to original directory.
                cd(d);
            end
            
            %% Generate object file of the C++/MATLAB bridge.
            brdIncPath = util.cppbridge.CppBridge.compilebridge(targetDir);
            
            %% Setup compile paramter.
            % CUDA path is hard coded for SDK v7.5, CUDA_INC_PATH is not
            % created for later CUDA SDKs, so we have to generate it by
            % ourselves.
            cudaIncPath = fullfile(getenv('CUDA_PATH_V7_5'), 'include');
            if exist(cudaIncPath, 'dir') ~= 7
                error(util.cppbridge.CppBridge.MSGID, ...
                      'CUDA_PATH is not properly configured.');
            else
                cudaIncPath = ['-I"', cudaIncPath, '"'];
            end
            
            %% Compile all the files.
            mexcuda('-v', ...
                    '-largeArrayDims', ...
                    '-c', ...
                    cudaIncPath, brdIncPath, ...
                    '*.cpp', '*.cu');
                
            %% Link all the files.
            mexcuda('-v', ...
                    '-output', mexName, ...
                    '*.obj');
        end

        function incPath = compilebridge(outDir)
            %COMPILEBRIDGE Compile the object file for C++/MATLAB bridge.
            %   P = COMPILEBRIDGE compiles the bridge object file and
            %   return the proper include path as P.

            %% Change directory into target folder.
            [brdDir, ~, ~] = fileparts(mfilename('fullpath'));
            localDir = cd(brdDir);

            %% Switch back the directory when COMPILEMEX is exited.
            cdLocalObj = onCleanup(@() cd(localDir));

            %% Generate the flag for include path.
            incPath = fullfile(brdDir, 'include');
            incPath = ['-I"', incPath, '"'];
            
            %% Compile the file.
            mex('-v', '-largeArrayDims', '-c', incPath, '*.cpp');
            
            %% Move the object file and calcaulte the include path.
            movefile('*.obj', outDir);
        end
    end
    
    methods (Abstract, Static, Access = protected)
        hdl = getfunchandle();
        %GETFUNCHANDLE Return the proper function handle.
        %   This function should be implemented by the childeren that
        %   intends to be instantiated. Since the MEX file may reside in a
        %   private folder, somewhere the parent function won't have direct
        %   access to.
    end
    
end
