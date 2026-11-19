function f = dda_figure(position, varargin)
% DDA_FIGURE  Figures for the DDA plots: tabs in one figure window, or separate windows.
%
%   f = dda_figure(position, 'Name', ...)   create a figure (name-value pairs as for figure)
%   dda_figure('docking', tf)               choose for the following figures:
%       true   tabs in MATLAB's common figure window (docked; the default)
%       false  separate windows of the given position/size
%
% The plot functions call dda_figure('docking', Input.plotting.docked) first.
% Docking also sets MATLAB's default window style, so figures made by helper
% functions (Q fit, FFT plots) become tabs too. In R2025a and later, giving a
% figure a Position undocks it, which is why the DDA plots create sized
% figures through this function.
if ischar(position) || isstring(position)
    if strcmpi(position, 'docking')
        docked = logical(varargin{1});
        setappdata(groot, 'DDA_figures_docked', docked);
        if docked, set(groot, 'DefaultFigureWindowStyle', 'docked');
        else,      set(groot, 'DefaultFigureWindowStyle', 'normal');
        end
        f = [];
        return
    end
    error('dda_figure:usage', 'Use dda_figure(position, ...) or dda_figure(''docking'', tf).');
end
docked = true;
if isappdata(groot, 'DDA_figures_docked'), docked = getappdata(groot, 'DDA_figures_docked'); end
if docked
    f = figure(varargin{:}, 'WindowStyle', 'docked');
else
    f = figure(varargin{:}, 'WindowStyle', 'normal', 'Position', position);
end
end
