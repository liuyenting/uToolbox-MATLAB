classdef Demo < handle
    %DEMO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access=private, SetObservable, AbortSet)
        value
    end
    
    methods
        function this = Demo(value)
            this.value = value;
            
            addlistener(this, 'value', 'PostSet', @gui.test.Demo.updateValue);
        end
        
        function this = setValue(this, value)
            this.value = value;
            fprintf('value is modified by setValue()\n');
        end
        
        function printValue(this)
            fprintf('value = %d\n', this.value);
        end
    end
    
    methods (Static, Access=private)
        function updateValue(source, event)
            this = event.AffectedObject;
            fprintf('value is updated to %d\n', this.value);
        end
    end
    
end

