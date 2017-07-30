classdef Model < handle
    %MODEL Represents the volumetric data.
    %
    %   TBA
    
    % observable properties, listeners are notified on change
    properties (SetObservable=true)
        voxelSize   % Voxel size along the X, Y and Z dimension.
        data        % Raw data.
        cursorPos   % Current cursor position in the volume.
    end
    
    % read-only properties
    properties (SetAccess=private)
        volumeSize  % Dimension of the volume.
    end
    
    % computable dependent properties
    properties (Dependent=true, SetAccess=private)
        ViewPlans   % Contains the XY, YZ and XZ view in order.
    end
    
    methods
        function this = Model()
            %CONSTRUCTOR
        end
        
        function plans = get.ViewPlans(this)
        end
    end
    
    % business logic
    
end

