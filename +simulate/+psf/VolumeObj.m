classdef VolumeObj < handle
    %VOLUMEOBJ Storing a 3D object, layer by layer.
    
    properties (Access = private)
        data;   % The 3D array that will hold the voxels. 
    end
    
    properties (GetAccess = private, SetAccess = immutable)
        nx;     % Size of X axis.
        ny;     % Size of Y axis.
        nz;     % Size of Z axis.
    end
    
    methods (Access = public)
        function obj = VolumeObj(varargin)
            % Initialize a zero 3D array from input using pixel unit.
            %   OBJ = VOLUME(DIM) creates the array by a 1-by-3 array in the
            %   order of X, Y and Z axis.
            %
            %   OBJ = VOLUME(W, H) creates the array with square planes of
            %   width W, and Z axis length of H.
            %
            %   OBJ = VOLUME(X, Y, Z) creates the array according to the
            %   definition X, Y and Z.
            
            switch nargin
                case 1
                    dimension = varargin{1};
                case 2
                    dimension = [varargin{1}, varargin{1}, varargin{2}];
                case 3
                    dimension = [varargin{1}, varargin{2}, varargin{3}];
                otherwise
                    error('PSFGenerator:Volume', 'Unknown definition for volume size.');
            end
            obj.data = zeros(dimension);
            [obj.nx, obj.ny, obj.nz] = size(obj.data);
        end
        
        function data = getdata(obj)
            data = obj.data;
        end
        
        function obj = setplane(obj, z, data)
            % Set the specified Z layer to DATA.
            
            obj.checkzbound('PSFGenerator:Volume:setplane', z);
            obj.checkdatadim('PSFGenerator:Volume:setplane', data);
            obj.data(:, :, z) = data;
        end
        
        function data = getplane(obj, z)
            % Get specified Z layer.
            
            obj.checkzbound('PSFGenerator:Volume:getplane', z);
            data = obj.data(:, :, z);
        end
    end
    
    methods (Access = private)
        function checkzbound(obj, msgID, z)
            if (z < 1) || (z > obj.nz)
                error(msgID, 'Z index exceeds boundary.');
            end
        end
        
        function checkdatadim(obj, msgID, data)
            if ~ismatrix(data)
                error(msgID, 'Data is not 2D.');
            elseif ~isequal([obj.nx, obj.ny], size(data))
                error(msgID, 'Dimension mismatched.');
            end
        end
    end
    
end

