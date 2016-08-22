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
    end

    methods (Static, Access = protected)
        function hdl = getfunchandle()
            hdl = str2func(local.Localize.MEX_NAME);
        end
    end
    
end
