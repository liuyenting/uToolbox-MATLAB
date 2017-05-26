classdef VolView < handle
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

        % hMultiView is a (2+k)x3 struct array, each row represents a handle for
        % different purpose, while each column represents XY/YZ/XZ multiview.
        % There are at least 2 types of axes - Raw and Crosshair.
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
            this.hMultiView = gobjects(2, 3);
            for i = 1:numel(this.hMultiView)
                this.hMultiView(i) = axes( ...
                    'Units', 'pixels', ...
                    'Position', [0, 0, 0, 0] ...
                );
            end

            % set layout properties
            this.fillRatio = 0.7;
            this.viewGap = 10;
            this.edgeGap = 40;

            % inject the data
            this.voxelSize = p.Results.VoxelSize;
            % default cursor position to the origin
            this.cursorPos = [200, 400, 80];

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
            % number of layers
            nl = size(this.hMultiView, 1);
            % iterate through the layers
            for l = 1:nl
                this.hMultiView(l, 1).Position = [e, e+Z+v, X, Y];
                this.hMultiView(l, 2).Position = [e+X+v, e+Z+v, Z, Y];
                this.hMultiView(l, 3).Position = [e, e, X, Z];
            end

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
                % iterate through the layers
                for l = 1:nl
                    % select the axes
                    axes(this.hMultiView(l, d));

                    % set the aspect ratio
                    set(gca, 'DataAspectRatioMode', 'manual');
                    set(gca, 'DataAspectRatio', [vsz(1:2), 1]);
                end

                % permute the voxel size for next setup
                %    XY             -> [px, py, 1] (px, py, pz)
                %    YZ             -> [py, pz, 1] (py, pz, px)
                %    XZ (tranposed) -> [pz, px, 1] (pz, px, py)
                vsz = permute(vsz, [2, 3, 1]);
            end
        end

        function this = updateAxes(this)
            disp('updateAxes()');

            % number of layers
            nl = size(this.hMultiView, 1);
            % iterate through the layers
            for l = 1:nl
                % XY
                axes(this.hMultiView(l, 1));
                    % X
                    xlabel('X', 'FontSize', 14);
                    set(gca, 'XAxisLocation', 'top');
                    set(gca, 'XTickLabel', []);
                    % Y
                    ylabel('Y', 'FontSize', 14);
                    set(gca, 'YTickLabel', []);

                % YZ
                axes(this.hMultiView(l, 2));
                    % X
                    xlabel('Z', 'FontSize', 14);
                    set(gca, 'XAxisLocation', 'top');
                    set(gca, 'XTickLabel', []);
                    % Y
                    set(gca, 'YTickLabel', []);

                % XZ
                axes(this.hMultiView(l, 3));
                    % X
                    set(gca, 'XTickLabel', []);
                    % Y
                    ylabel('Z', 'FontSize', 14);
                    set(gca, 'YTickLabel', []);
            end
        end

        function this = updateMultiView(this)
            disp('updateMultiView()');

            this.updateRawData();
            this.updateCrosshair();
            this.updateAdditionalLayers();
        end

        function this = updateRawData(this)
            disp('updateRawData()');

            plotter = @imagesc;

            for d = 1:3
                axes(this.hMultiView(1, d));
                % target dimension
                td = (3-d)+1;
                A = util.ndslice(this.data, td, this.cursorPos(td));
                % transpose if necessary
                if d == 3
                    A = A.';
                end
                plotter( ...
                    A, ...
                    'HitTest', 'off' ... % use underlying axes for events
                );
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
                axes(this.hMultiView(2, d));

                % remove the background
                %set(gca, 'Color', 'none');

                % draw the crosshair
                hold on;
                % vertical
                x = [c(1),  c(1)];
                y = [   0, sz(2)];
                line(x, y, 'Color', 'yellow');
                % horizontal
                x = [   0, sz(1)];
                y = [c(2),  c(2)];
                line(x, y, 'Color', 'yellow');
                hold off;

                sz = permute(sz, [2, 3, 1]);
                c = permute(c, [2, 3, 1]);
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
    end
end
