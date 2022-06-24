function [G_N, N_A, N] = GrainSize_E112_SaltikovPlanimetric(ebsd, varargin)
% Saltikov's Planimetric: Count of grains in a test rectangle
% Grains completely enclosed count as 1, those intercepted by the rectangle
% count by half, and the corners count as one-quarter.

% make a rectangle
offset = 0.02; % 2pct inset from edges
xres = (max(ebsd.x) - min(ebsd.x)) / length(ebsd.x);
yres = (max(ebsd.y) - min(ebsd.y)) / length(ebsd.y);
xinset = numel(ebsd.x) * offset * xres;
yinset = numel(ebsd.y) * offset * yres;

polygon = [min(ebsd.x)+xinset, min(ebsd.y)+yinset;
           min(ebsd.x)+xinset, max(ebsd.y)-yinset;
           max(ebsd.x)-xinset, max(ebsd.y)-yinset;
           max(ebsd.x)-xinset, min(ebsd.y)+yinset];

[G_N, ~, N_A, N, ~, ~, inside_grains2, edge_grains2, ~] = grainsize_areas_planimetric(ebsd, polygon, varargin{:});

if ismember('PlotResults',varargin)
    plot(ebsd, ebsd.orientations); hold on
    plot(edge_grains2.boundary, 'linewidth', 2, 'lineColor', 'black');
    plot(inside_grains2.boundary, 'linewidth', 3, 'lineColor', 'white');
    plot(polygon(:,1), polygon(:,2), 'k', 'linewidth', 3)
end % parse varargin

end