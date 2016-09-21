classdef ProgressBar
    %PROGRESSBAR Generate and manipulate a progress bar.
    %   Detailed explanation goes here
    
    properties (Constant, Access = private, Hidden = true)
        MSGID = 'gui:progressbar';
    end
    
    properties (Constant, Access = private)
        REG_WIDTH = 512;    % Width of the bar region.
        REG_HEIGHT = 64;    % Height of a single bar region.
        BAR_WIDTH = 384;    % Width of the bar.
        BAR_HEIGHT = 24;    % Height of the bar.
    end
    
    methods (Static)
        function h = create(title)
            %CREATE Create a new progress bar.
            %   H = CREATE() creates a progress bar object.
            %   H = CREATE(TITLE) creates an object with figure name set to
            %   TITLE.
            
            %  Create and then hide the UI as it is being constructed.
            if nargin == 0
                h = figure('Visible', 'off');
            else
                if nargin > 1
                    warning(gui.ProgressBar.MSGID, ...
                            'Additional arguments ignored.');
                end
                h = figure('Visible', 'off', ...
                           'Name', title, ...
                           'NumberTitle', 'off');
            end
            
            h.WindowStyle = 'normal';
            
            % Hide unused features.
            h.DockControls = 'off';
            h.MenuBar = 'none';
            h.ToolBar = 'none';
            h.Resize = 'off';
            
            % Set default internal data.
            data.NumberOfBars = 0;
            guidata(h, data);
            
            % Configure the default size.
            h.OuterPosition = gui.ProgressBar.calcposvec(0);
            
            h.Visible = 'on';
        end
        
        function h = addBar(h, nbar, varargin)
            %ADDBAR Add bars to the progress bar figure.
            
            if nbar < 0
                error(gui.ProgressBar.MSGID, ...
                      'Only positive number of bars is valid.');
            elseif nbar == 0
                return;
            end
            
            h.Visible = 'off';
            
            for i = 1:nbar
                desc = varargin{2*i-1};
                val = varargin{2*i};
                
                if isempty(desc)
                    % Place the progress bar in the center.
                end
            end
            
            % Update the counter.
            data = guidata(h);
            data.NumberOfBars = data.NumberOfBars + nbar;
            guidata(h, data);
            
            % Resize the window.
            h.OuterPosition = gui.ProgressBar.calcposvec(nbar);
            
            h.Visible = 'on';
        end
        
        function h = setText(h, ibar, desc)
        end
        
        function h = setProgress(h, ibar, val)
        end
    end
    
    methods (Static, Access = private)
        function v = calcposvec(nbar)
            %CALCPOSVEC Calculate the screen position vector.
            
            % Default to one if there is no progress bar.
            if nbar < 0
                error(gui.ProgressBar.MSGID, ...
                      'Number of progress bar should be positive.');
            elseif nbar == 0
                nbar = 1;
            end
            
            fw = gui.ProgressBar.REG_WIDTH;
            fh = gui.ProgressBar.REG_HEIGHT * nbar;
            
            [sw, sh] = util.screensize();
            fl = (sw-fw)/2;
            fb = (sh-fh)/2;
            
            v = [fl, fb, fw, fh];
        end
    end
end

