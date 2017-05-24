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

    properties (Access=private, Hidden=true)
        hFigure;        % Handle to the actual figure.
        hAxes;          % Handles to XY/YZ and XZ axes.
        hPreview;       % Handles to the 3-D preview axes.
        hGraphics;      % Handles to the plotted contents.
    end

    properties (SetAccess=protected, GetAccess=protected, Hidden=true)
        axSep;          % Interval between axes.
        edgeSep;        % Padding length between the axes and window.
        fillRatio;      % Filling ratio for the figure versus dispaly.

        voxSize;        % Voxel size along X, Y and Z dimension.
        data;           % Volume to display.
        volSize;        % Dimension of the volume.

        cursorPos;      % Current cursor position in the data.
        %ROI
    end

    methods
        function this = VolView(varargin)
            p = inputParser;
            % only 3-D data is allowed
            addOptional(p, 'Data', [], @(x)(~isempty(x) && (ndims(x)==3)));
            % voxels are default to be isotropic
            addParameter(p, 'VoxelSize', [1, 1, 2.5], @(x)(isnumeric(x)));
            % use variable name as the default title
            addParameter(p, 'Title', inputname(1), @(x)(ischar(x)));
            parse(p, varargin{:});

            % populate the figure
            this.hFigure = figure( ...
                'Name', p.Results.Title, ...
                'NumberTitle', 'off', ...
                'Visible', 'off' ...
            );

            % construct the components
            this.hAxes = gobjects(1, 3);
            for i = 1:3
                this.hAxes(i) = axes( ...
                    'Units', 'pixels', ...
                    'Position', [0, 0, 0, 0] ...
                );
            end
            this.hPreview = gobjects;
            this.hGraphics = gobjects(1, 3);

            % set figure related variables
            this.axSep = 10;
            this.edgeSep = 40;
            this.fillRatio = 0.7;

            this.voxSize = p.Results.VoxelSize;
            data = p.Results.Data;
            if ~isempty(data)
                this.show(data);
            end
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

            this.data = data;

            % assess the dimension
            sz = size(data);
            %TODO: process data with multi-timepoints and channels
            this.volSize = sz;

            if gui.volview.VolView.DEBUG
                fprintf('volume = %dx%dx%d\n', sz(1), sz(2), sz(3));
            end

            this.updateAxesSize();
            this.updateAxesContent();
            this.updateAxesAnnotation();
        end

        function this = updateAxesSize(this)
            this.hFigure.Visible = 'off';

            % retrieve volume size, compensated by voxel size to isotropic
            sz = this.volSize .* this.voxSize;
            % set window size, [X+Z+axSep+2*edgeSep, Y+Z+axSep+2*edgeSep]
            xyzSz = sz(1:2) + sz(3);

            % expand the window to fill the window
            scSz = util.screensize;
            % only fill to designated percentage by expanding the axes
            scSz = scSz * this.fillRatio;
            ratio = scSz./xyzSz;
            % use the smaller ratio to expand
            ratio = min(ratio);
            xyzSz = xyzSz * ratio;

            % filler size
            filler = this.axSep + 2*this.edgeSep;
            % update the size
            winSz = xyzSz + filler;
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
            X = sz(1)*ratio;
            Y = sz(2)*ratio;
            Z = sz(3)*ratio;
            this.hAxes(1).Position = [e, e+Z+a, X, Y];
            this.hAxes(2).Position = [e+X+a, e+Z+a, Z, Y];
            this.hAxes(3).Position = [e, e, X, Z];

            this.hFigure.Visible = 'on';
        end

        function this = updateAxesContent(this, pos)
            %UPDATEAXESCONTENT Updates the view.
            %
            %   THIS = UPDATEAXESCONTENT(THIS) set the slices to center of
            %   the volume.
            %   THIS = UPDATEAXESCONTENT(THIS, POS) set the slices to the
            %   specified location and update the internal variable.

            if nargin == 1
                pos = this.volSize/2;
                pos = floor(pos);
            end

            plotter = @imagesc;
            % XY
            axes(this.hAxes(1));
            this.hGraphics(1) = plotter(this.data(:, :, pos(3)));

            % YZ
            axes(this.hAxes(2));
            this.hGraphics(2) = plotter(squeeze(this.data(:, pos(2), :)));

            % XZ
            axes(this.hAxes(3));
            this.hGraphics(3) = plotter(squeeze(this.data(pos(1), :, :)).');

            % apply colormap
            colormap(gray);

            this.cursorPos = pos;
        end

        function this = updateAxesAnnotation(this)
            px = this.voxSize(1);
            py = this.voxSize(2);
            pz = this.voxSize(3);

            % XY
            axes(this.hAxes(1));
            xlabel('X', 'FontSize', 14);
            ylabel('Y', 'FontSize', 14);
            set(gca, 'XAxisLocation', 'top');
            set(gca, 'XTickLabel', []);
            set(gca, 'YTickLabel', []);
            set(gca, 'DataAspectRatioMode', 'manual');
            set(gca, 'DataAspectRatio', [px, py, 1]);

            % YZ
            axes(this.hAxes(2));
            xlabel('Z', 'FontSize', 14);
            set(gca, 'XAxisLocation', 'top');
            set(gca, 'XTickLabel', []);
            set(gca, 'YTickLabel', []);
            set(gca, 'DataAspectRatioMode', 'manual');
            set(gca, 'DataAspectRatio', [py, pz, 1]);

            % XZ (transposed)
            axes(this.hAxes(3));
            ylabel('Z', 'FontSize', 14);
            set(gca, 'XTickLabel', []);
            set(gca, 'YTickLabel', []);
            set(gca, 'DataAspectRatioMode', 'manual');
            set(gca, 'DataAspectRatio', [pz, px, 1]);
        end
    end

    methods (Static, Access=private, Hidden=true)

    end

end
