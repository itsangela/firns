classdef datefig < handle
    % DATEFIG A date figure handle class object.
    %   OBJ = DATEFIG creates a date figure object.
    %
    %   [H,OBJ,DATA] = DATEFIG.PLOT(...) creates a date figure object OBJ
    %   and plots DATA with date ticks on the x-axis.
    %   DATEFIG.PLOT(OBJ,...) creates a plot in OBJ, which must be a
    %   datefig class object.
    %
    %   DATEFIG.SCATTER(...) creates a date figure object OBJ with a
    %   scatter plot and places dateticks on the x-axis of a scatter plot.
    %   DATEFIG.SCATTER(OBJ,...) creates scatter plot in OBJ which must be
    %   a datefig class object.
    %
    %   [AX,H1,H2,OBJ,DATA] = DATEFIG.PLOTYY(X1,Y1,X2,Y2) creates a date
    %   figure object OBJ with PLOTYY and places dateticks on the x-axis of
    %   a scatter plot.
    %   DATEFIG.PLOTYY(OBJ,...) creates PLOTYY in OBJ which must be a
    %   datefig class object.
    %
    % See also DATETICK, ZOOM, PAN, PLOT/PLOT, SCATTER/SCATTER,
    %   PLOTYY/PLOTYY, FIGURE, LINKAXES
    %
    % Reference:
    % http://www.mathworks.com/help/releases/R2011a/techdoc/ref/figure_props.html#ResizeFcn

    % Mark Mikofski, SunPower Corp.
    % Version 1-0, 2011-07
    % Version 1-1, 2011-07-28
    % 1. fix error, plot was only accepting pairs of inputs, could not
    %   input string for color/shape, single argument, etc.
    % 2. add overloaded method for scatter
    % Version 1-2, 2011-07-28
    % add overloaded method for plotyy, fix error that cause uitool to
    % throw error on datetick since eventData.Axes is vector and datetick
    % only accepts scalar, fix error with reset zoom/pan only rescales one
    % of axes
    
    properties (Constant)
        version = 1.2 % version of code
    end
    properties (SetAccess = immutable)
        figureHandle = [] % figure handle
        zoomHandle = [] % zoom handle
        panHandle = [] % pan handle
    end
    methods
        function df = datefig
            df.figureHandle = figure;
            set(df.figureHandle,'ResizeFcn',@datefig.myResizeFcn)
            df.zoomHandle = zoom(df.figureHandle);
            df.panHandle = pan(df.figureHandle);
            set(df.zoomHandle,'ActionPostCallback',@datefig.myZoomPanFcn);
            set(df.panHandle,'ActionPostCallback',@datefig.myZoomPanFcn);
        end % constructor
    end
    methods (Static)
        function [h,obj,data] = plot(varargin)
            % DATEFIG.PLOT Plot in a date figure.
            %   [H,OBJ,DATA] = DATEFIG.PLOT(...) creates a date figure
            %   object OBJ and plots DATA with date ticks on the x-axis.
            %   DATEFIG.PLOT(OBJ,...) creates a plot in OBJ, which must be
            %   a datefig class object.
            if isa(varargin{1},'datefig') % checks if first argument is a DATEFIG object
                obj = varargin{1};
                figure(obj.figureHandle); % makes OBJ the current figure
                data = varargin(2:end);   % data to pass to PLOT
            else
                obj = datefig;
                data = varargin;
            end
            h = plot(data{:});
            datetick('x')
            xlabel('Time')
            grid
        end
        function [h,obj,data] = scatter(varargin)
            % DATEFIG.SCATTER Create a scatter plot in a date figure.
            %   [H,OBJ,DATA] = DATEFIG.SCATTER(...) creates a date figure
            %   object OBJ with a scatter plot and places dateticks on the
            %   x-axis of a scatter plot.
            %   DATEFIG.SCATTER(OBJ,...) creates scatter plot in OBJ which
            %   must be a datefig class object.
            if isa(varargin{1},'datefig') % checks if first argument is a DATEFIG object
                obj = varargin{1};
                figure(obj.figureHandle); % makes OBJ the current figure
                data = varargin(2:end);   % data to pass to PLOT
            else
                obj = datefig;
                data = varargin;
            end
            h = scatter(data{:});
            datetick('x')
            xlabel('Time')
            grid
        end
        function [ax,h1,h2,obj,data] = plotyy(varargin)
            % DATEFIG.PLOTYY Create a 2-axes plot in a date figure.
            %   [AX,H1,H2,OBJ,DATA] = DATEFIG.PLOTYY(X1,Y1,X2,Y2) creates a
            %   date figure object OBJ with PLOTYY and places dateticks on
            %   the x-axis of a scatter plot.
            %   DATEFIG.PLOTYY(OBJ,...) creates PLOTYY in OBJ which must be
            %   a datefig class object.
            if isa(varargin{1},'datefig') % checks if first argument is a DATEFIG object
                obj = varargin{1};
                figure(obj.figureHandle); % makes OBJ the current figure
                data = varargin(2:end);   % data to pass to PLOT
            else
                obj = datefig;
                data = varargin;
            end
            [ax,h1,h2] = plotyy(data{:});
            linkaxes(ax,'x')
            datetick(ax(1),'x','keeplimits'),datetick(ax(2),'x','keeplimits')
            xlabel('Time')
            grid
        end
        function myZoomPanFcn (figureHandle,eventData) %#ok<INUSD>
            % zoomHandle.ActionPostCallback requires figureHandle arg.
            subPlots = allchild(figureHandle);
            subPlots = subPlots(strcmp(get(subPlots,'Type'),'axes'));
            if ~isempty(subPlots)
                subPlots = subPlots(~strcmp(get(subPlots,'Tag'),'legend'));
            end
            nSubPlots = length(subPlots);
            for n = 1:nSubPlots
                datetick(subPlots(n),'x','keeplimits');
            end
            % Don't use eventData.Axes becuase:
            % 1. For zoom/pan on PLOTYY, eventData.Axes is a vector, and
            % DATETICK will fail since it does not accept vector of axes.
            % 2. For reset zoom/pan on PLOTYY, eventData.Axes only contains
            % one axes, and so not all axes will be rescaled.
        end
        function myResizeFcn(figureHandle,eventData) %#ok<INUSD>
            % figureHandle.ResizeFcn requires eventData arg.
            subPlots = allchild(figureHandle);
            subPlots = subPlots(strcmp(get(subPlots,'Type'),'axes'));
            if ~isempty(subPlots)
                subPlots = subPlots(~strcmp(get(subPlots,'Tag'),'legend'));
            end
            nSubPlots = length(subPlots);
            for n = 1:nSubPlots
                datetick(subPlots(n),'x','keeplimits');
            end
        end
    end
    
end