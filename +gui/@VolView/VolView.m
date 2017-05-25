classdef VolView < handle
    %VOLVIEW Display volumetric data set.
    %
    %   H = VOLVIEW() creates a default viewer.
    %   H = VOLVIEW(DATA) shows render the volume of the data in XY/YZ/XZ
    %   view.
    %   H = VOLVIEW(..., 'Title', NAME) set the title of the viewer instead
    %   of default to the variable name.

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
    end

    %% Constructor and destructor
    methods
        function this = VolView(varargin)
            %CONSTRUCTOR Create a template volume viewer object.

        end

        function delete(this)
            %DESTRUCTOR Free all the resources.

            % close remaining figure
            h = this.hFigure;
            if ~isempty(h) && ishandle(h)
                close(h);
            end
        end
    end

    %% Public functions
    methods

    end
end
