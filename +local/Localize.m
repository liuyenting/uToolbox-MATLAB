classdef Localize < util.cppbridge.CppBridge
    %LOCALIZE Perform sub-pixel localization.

    properties (Constant, Access = private, Hidden = true)
        MSG_ID = 'local:localize:matlab';
        MEX_NAME = 'localizemex';
    end

    properties (Access = private)
    end

    methods
        function this = Localize()
            this@util.cppbridge.CppBridge(mfilename('fullpath'), ...
                                          local.Localize.MEX_NAME);
        end
        
        function this = setparm(this, parm)
            %SETPARM Set the fitting parameters.
            
            % Make sure the fields exists.
            if ~isfield(parm, 'FittingSize')
                error(local.Localize.MSG_ID, ...
                      'Missing the size of fitting ROI.');
            elseif ~isfield(parm, 'ThresholdRatio')
                error(local.Localize.MSG_ID, ...
                      'Missing the ratio of filter threshold.');
            end
            
            if (parm.FittingSize < 0) || (parm.ThresholdRatio < 0)
                error(local.Localize.MSG_ID, ...
                      'Parameters should be positive only.');
            end
            
            % Coerce the data type.
            parm.FittingSize = uint16(parm.FittingSize);
            
            if parm.ThresholdRatio < 0.5
                warning(local.Localize.MSG_ID, ...
                        'Coerce the ratio to 0.5');
            elseif parm.ThresholdRatio > 2
                warning(local.Localize.MSG_ID, ...
                        'Coerce the ratio to 2.');
            end
            parm.ThresholdRatio = single(parm.ThresholdRatio);
            
            this.mexHandle('setparm', this.objectHandle, ...
                           parm);
        end
        
        function this = loadstack(this, im)
            %LOADSTACK Load the image stack to be processed.
           
            
        end
    end

    methods (Static, Access = protected)
        function hdl = getfunchandle()
            hdl = str2func(local.Localize.MEX_NAME);
        end
    end
    
end
