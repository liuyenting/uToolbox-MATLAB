classdef VolView < handle
    %VOLVIEW Display volumetric data set.
    %
    %   H = VOLVIEW() creates a default viewer.
    %   H = VOLVIEW(DATA) shows render the volume of the data in XY/YZ/XZ
    %   view.
    %   H = VOLVIEW(..., 'Title', NAME) set the title of the viewer instead
    %   of default to the variable name.
    %
    %   See also: TBA
    
    properties (Constant, Access=private, Hidden=true)
        MSGID = 'gui:volview';
        DEBUG = false;
    end 
    
    properties (SetAccess=protected, GetAccess=protected, Hidden=true)
        hFigure;      % Handle to the actual figure.
        hAxes;       % Handles to XY/YZ and XZ axes.
        data;           % Volume to display. 
        volSize;        % Dimension of the volume.
        
        axSep;          % Interval between axes.
        edgeSep;        % Padding length between the axes and window.
    end
    
    methods
        function this = VolView(varargin)
            p = inputParser;
            % only 3-D data is allowed
            addOptional(p, 'Data', [], @(x)(~isempty(x) && (ndims(x)==3)));
            % use variable name as the default title
            addParameter(p, 'Title', inputname(1), @(x)(ischar(x)));
            parse(p, varargin{:});
            
            % populate the figure
            this.hFigure = figure( ...
                'Name', p.Results.Title, ...
                'NumberTitle', 'off', ...
                'Visible', 'off' ...
            );
            
            this.setData(p.Results.Data);
            
            % construct the components
            this.hAxes = gobjects(1, 3);
            for i = 1:3
                this.hAxes(i) = axes( ...
                    'Units', 'pixels', ...
                    'Position', [0, 0, 0, 0] ...
                );
            end
            
            % set figure related variables
            this.axSep = 20;
            this.edgeSep = 40;
        end
        
        function this = show(this, data)
            if isempty(data)
                error(gui.volview.VolView.MSG_ID, 'Input must be valid!');
            else
                this.setData(data);
            end
        end
        
        function delete(this)
            h = this.hFigure;
            % close the figure if exists
            if ~isempty(h) && ishandle(h)
                close(h);
            end
        end
    end
    
    methods (Access=private, Hidden=true)
        function this = setData(this, data)
            %SETDATA Set the data to display.
            
            if ~isempty(data)
                this.data = data;
                
                % assess the dimension
                sz = size(data);
                %TODO: process data with multi-timepoints and channels
                this.volSize = sz;
                
                if gui.volview.VolView.DEBUG
                    fprintf('volume = %dx%dx%d\n', sz(1), sz(2), sz(3));
                end
                
                this.updateAxesSize();
            end
        end
        
        function this = updateAxesSize(this)
            this.hFigure.Visible = 'off';
            
            % retrieve volume size
            sz = this.volSize;
            % set window size, [X+Z+axSep+2*edgeSep, Y+Z+axSep+2*edgeSep]
            winSz = sz(1:2) + sz(3);
            winSz = winSz + this.axSep + 2*this.edgeSep;
            % update the size
            p = this.hFigure.Position;
            this.hFigure.Position = [p(1:2), winSz];
            % center the window
            movegui(this.hFigure ,'center');
            
            % +----+----+
            % | XY | YZ |
            % +----+----+
            % | XZ |    |
            % +----+----+
            % 
            %   XY
            %       Left = e
            %       Bottom = e + Z + a
            %       Width = X
            %       Height = Y
            %   YZ
            %       Left = e + X + a
            %       Bottom = e + Z + a
            %       Width = Z
            %       Height = Y
            %   XZ
            %       Left = e
            %       Bottom = e
            %       Width = X
            %       Height = Z
            e = this.edgeSep;
            a = this.axSep;
            X = sz(1);
            Y = sz(2);
            Z = sz(3);
            this.hAxes(1).Position = [e, e+Z+a, X, Y];
            this.hAxes(2).Position = [e+X+a, e+Z+a, Z, Y];
            this.hAxes(3).Position = [e, e, X, Z];
            
            this.hFigure.Visible = 'on';
        end
    end
    
end

