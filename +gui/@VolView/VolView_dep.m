classdef VolView_dep < handle
    %VOLVIEW Display volumetric data set.
    %
    %   H = VOLVIEW() creates a default viewer.
    %   H = VOLVIEW(DATA) shows render the volume of the data in XY/YZ/XZ
    %   view.
    %   H = VOLVIEW(..., PARAM) allows detail controls of the internal
    %   parameters.
    %
    %   Parameters
    %   ----------
    %   'VoxelSize'     Set the voxel size if the data is anisotropic.
    %   'Title'         Title of the screen, default to DATA variable name.

    %% Book-keeping variables
    properties (Access=private, Hidden=true)
        % hFigure holds the handle to volume viwer's root figure object.
        hFigure;

        % hMultiView is a Nx3 graphic object array, each column represents
        % XY, YZ and XZ view respectively. Default row is defined as
        % 	1) Axes
        %   2) Raw data
        %   3..) Additional data
        hMultiView;

        % hPreview holds the handles for the overview of current position in the
        % volumetric data.
        hPreview;

        % hListener holds an array that contains all the event listener applied
        % in the constructor. Used by the destructor to unregister all the user
        % events.
        hListener;
    end

    %% Layout configurations
    properties (SetAccess=protected, GetAccess=public)
        fillRatio;      % Ratio of the entire viewer respective to the screen.
        viewGap;        % Gaps (px) between the views.
        edgeGap;        % Gaps (px) between the views and the edges.
    end

    %% Data
    properties (SetAccess=protected, GetAccess=public, SetObservable, AbortSet)
        voxelSize;      % Voxel size along the X, Y and Z dimension.
        volumeSize;     % Dimension of the volume.
        %TODO: attach volume size variation to axes poisition update function.

        data;           % Raw data.
        %TODO: attach data change to hGraphics update callback

        cursorPos;      % Current cursor position in the
        %TODO: update preview and boundary check by cursorPos change
    end

    %% Constructor and destructor
    methods
        function this = VolView(varargin)
            %CONSTRUCTOR Create a template volume viewer object.

            p = inputParser;
            % only 3-D data is allowed
            addOptional(p, 'Data', [], @(x)(~isempty(x) && (ndims(x)==3)));
            % voxels are default to be isotropic
            addParameter(p, 'VoxelSize', [1, 1, 1], @(x)(isnumeric(x)));
            % use variable name as the default title
            addParameter(p, 'Title', inputname(1), @(x)(ischar(x)));
            parse(p, varargin{:});

            % generate the figure
            this.hFigure = figure( ...
                'Name', p.Results.Title, ...
                'NumberTitle', 'off', ...
                'Visible', 'off' ...
            );

            % populate the axes for raw data and crosshairs
            this.hMultiView = gobjects(3, 3);
            for d = 1:3
                this.hMultiView(1, d) = axes( ...
                    'Units', 'pixels', ...
                    'Position', [0, 0, 0, 0], ...
                    'YDir', 'reverse', ...  % default to image style
                    'NextPlot', 'add' ...   % keep parent values
                );
            end
            % add the dimension tag
            this.hMultiView(1, 1).Tag = 'XY';
            this.hMultiView(1, 2).Tag = 'YZ';
            this.hMultiView(1, 3).Tag = 'XZ';

            % set layout properties
            this.fillRatio = 0.7;
            this.viewGap = 10;
            this.edgeGap = 40;

            % inject the data
            this.voxelSize = p.Results.VoxelSize;
            % default cursor position to the origin
            this.cursorPos = [900, 600, 65];

            % attach the listener
            propName = {'voxelSize', 'volumeSize', 'data', 'cursorPos'};
            np = numel(propName);
            for i = 1:np
                lh = addlistener( ...
                    this, propName{i}, ...
                    'PostSet', @gui.VolView.propertyChangeEvents ...
                );
                this.hListener = [this.hListener, lh];
            end
            this.hFigure.WindowButtonDownFcn = @gui.VolView.mouseDown;

            % inject the data
            this.data = p.Results.Data;
        end

        function delete(this)
            %DESTRUCTOR Free all the resources.

            % remove registered user event handlers
            lh = this.hListener;
            if ~isempty(lh)
                for l = lh
                    delete(l);
                end
            end

            % close remaining figure
            fh = this.hFigure;
            if ~isempty(fh) && ishandle(fh)
                close(fh);
            end
        end
    end

    %% Public functions
    methods
        this = show(this, data)
        this = setCursor(this, pos)
    end
    
    %% Setters and Getters
    methods
        function set.volumeSize(this, sz)
            % Note: MATLAB layouts 3-D matrix as (Y, X, Z)
            sz([1, 2]) = sz([2, 1]);
            this.volumeSize = sz;
        end
    end
    
    %% Primary callback function for events
    methods (Static, Access=private)
        function propertyChangeEvents(source, event)
            % primary object
            this = event.AffectedObject;

            switch source.Name
                case 'voxelSize'
                    disp('update "voxelSize"');

                    this.updateLayout();
                    this.updateAspectRatio();
                case 'volumeSize'
                    disp('update "volumeSize"');

                    %DEBUG...
                    sz = this.volumeSize;
                    fprintf('volume = %dx%dx%d\n', sz(1), sz(2), sz(3));
                    %...DEBUG
                    this.updateLayout();
                    this.updateAxes();
                case 'data'
                    disp('update "data"');

                    % update volume size
                    this.volumeSize = size(this.data);
                    % update the display
                    this.updateMultiView();
                    this.updateAxes();
                case 'cursorPos'
                    disp('update "cursorPos"');

                    this.updateCrosshair();
                    this.updateMultiView();
                    this.updateAxes();
            end
        end
        
        function mouseDown(source, event)
            % primary object
            this = event.AffectedObject;
            
            h = gca;
            cNew = get(h, 'CurrentPoint');
            disp(['clicked @ ', h.Tag, '!']);
            
            %cOld = this.
            % generate new cursor position
            switch h.Tag
                case 'XY'
                    
                case 'YZ'
                case 'XZ'
            end
        end
    end

    %% Private functions
    methods (Access=private)
        function this = updateLayout(this)
            disp('updateLayout()');

            this.hFigure.Visible = 'off';

            % retrieve volume size, compensated by voxel size to isotropic
            sz = this.volumeSize .* this.voxelSize;
            % set plot multi-view size
            xyzSz = sz(1:2) + sz(3);

            % filler size
            filler = this.viewGap + 2*this.edgeGap;
            % get the window size
            scSz = util.screensize;
            % only fill to designated percentage by expanding the axes
            scSz = scSz * this.fillRatio;
            % expand the multi-view to maximum plausible screen area...
            ratio = (scSz-filler) ./ xyzSz;
            % ... use the smaller ratio to expand
            ratio = min(ratio);
            xyzSz = xyzSz * ratio;

            % update the window size
            winSz = xyzSz + filler;
            p = this.hFigure.Position;
            this.hFigure.Position = [p(1:2), winSz];
            % center the window to screen
            movegui(this.hFigure ,'center');

            % +----+----+
            % | XY | YZ |
            % +----+----+
            % | XZ |    |
            % +----+----+
            %
            %   XY
            %       Left = e
            %       Bottom = e + Z + v
            %       Width = X
            %       Height = Y
            %   YZ
            %       Left = e + X + v
            %       Bottom = e + Z + v
            %       Width = Z
            %       Height = Y
            %   XZ
            %       Left = e
            %       Bottom = e
            %       Width = X
            %       Height = Z
            e = this.edgeGap;
            v = this.viewGap;
            X = sz(1)*ratio;
            Y = sz(2)*ratio;
            Z = sz(3)*ratio;
            
            % apply to the axes
            this.hMultiView(1, 1).Position = [e, e+Z+v, X, Y];
            this.hMultiView(1, 2).Position = [e+X+v, e+Z+v, Z, Y];
            this.hMultiView(1, 3).Position = [e, e, X, Z];

            this.hFigure.Visible = 'on';
        end

        function this = updateAspectRatio(this)
            disp('updateAspectRatio()');

            % number of layers
            nl = size(this.hMultiView, 1);
            % copied voxel size to avoid tampering
            vsz = this.voxelSize;

            % iterate through XY, YZ, XZ
            for d = 1:3
                % select the axes
                h = this.hMultiView(1, d);

                % set the aspect ratio
                h.DataAspectRatioMode = 'manual';
                h.DataAspectRatio = [vsz(1:2), 1];

                % permute the voxel size for next setup
                %    XY             -> [px, py, 1] (px, py, pz)
                %    YZ             -> [py, pz, 1] (py, pz, px)
                %    XZ (tranposed) -> [pz, px, 1] (pz, px, py)
                vsz = circshift(vsz, -1);
            end
        end

        function this = updateAxes(this)
            %UPDATEAXES Update underlying axes render configurations.
            %   
            %   TBA
            %   Remove ticks and set appropriate labels.
    
            % XY
            h = this.hMultiView(1, 1);
                % X
                xlabel(h, 'X', 'FontSize', 14);
                h.XAxisLocation = 'top';
                h.XTick = [];
                % Y
                ylabel(h, 'Y', 'FontSize', 14);
                h.YTick = [];

            % YZ
            h = this.hMultiView(1, 2);
                % X
                xlabel(h, 'Z', 'FontSize', 14);
                h.XAxisLocation = 'top';
                h.XTick = [];
                % Y
                h.YTick = [];

            % XZ
            h = this.hMultiView(1, 3);
                % X
                h.XTick = [];
                % Y
                ylabel(h, 'Z', 'FontSize', 14);
                h.YTick = [];
        end

        function this = updateMultiView(this)
            disp('updateMultiView()');

            this.updateRawData();
            this.updateCrosshair();
            this.updateAdditionalLayers();
        end

        function this = updateRawData(this)
            plotter = @imagesc;
            
            % MATLAB layouts 3-D array as..    (Y, X, Z)
            % (XY, YZ, XZ) slices defined by.. (Z, Y, X) 
            % Internal size definition is..    (X, Y, Z)
            %
            % So, the slices are defined by (in permuted index) [3, 2, 1],
            % but we have to swap the X/Y index, so they are called by [3,
            % 1, 2], the worst possible permutation.
            indSel = [3, 1, 2];
        
            % iterate through views
            for d = 1:3
                % target dimension
                td = (3-d)+1;
                % selected layer is controlled by pre-calculated indices
                ci = indSel(d);
                A = util.ndslice(this.data, td, this.cursorPos(ci));
                % transpose if necessary
                if d == 3
                    A = A.';
                end
                
                h = this.hMultiView(1, d);
                this.hMultiView(2, d) = plotter( ...
                    h, ...                  % target axes
                    A, ...                  % the sliced data
                    'HitTest', 'off' ...    % don't process events 
                );
                % tighten the axis
                axis(h, 'tight');
            end

            % apply colormap
            colormap(gray);
        end

        function this = updateCrosshair(this)
            % size of the volume
            sz = this.volumeSize;
            % cursor
            c = this.cursorPos;
            
            % iterate through the dimensions
            for d = 1:3
                h = this.hMultiView(1, d);

                % vertical
                x = [c(1),  c(1)];
                y = [   0, sz(2)];
                if d ~= 1
                    t = x; x = y; y = t;
                end
                line( ...
                    h, ...
                    x, y, ...
                    'Color', 'yellow', ...
                    'LineWidth', 1, ...
                    'HitTest', 'off' ...
                );
                
                % horizontal
                x = [   0, sz(1)];
                y = [c(2),  c(2)];
                if d ~= 1
                    t = x; x = y; y = t;
                end
                line( ...
                    h, ...                  % select view axes
                    x, y, ...               % start/end matrix
                    'Color', 'yellow', ...  % line apperances
                    'LineWidth', 1, ...
                    'HitTest', 'off' ...    % don't process events
                );
                
                sz = circshift(sz, -1);
                c = circshift(c, -1);
            end
        end

        function this = updateAdditionalLayers(this)
            % number of layers
            nl = size(this.hMultiView, 1);

            % skip this process if only default layer exists
            if nl == 2
                return;
            end

            %TODO
        end
        
        function this = updateCursorPos(this)
            
        end
    end
end
