function [G_N, N_A, N] = GrainSize_E112_JeffriesPlanimetric(ebsd, varargin)
% Jeffries' Planimetric: Count of grains in a test circle
% Grains completely enclosed count as 1, those intercepted by the circle count by half. 

% approximate a circle as a polygon
offset = 0.02; % 2pct inset from edges
xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x));
ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y));
thetas = 0.0:pi/100.0:2.0*pi;
xres = 2.0 * xcenter / length(ebsd.x);
yres = 2.0 * ycenter / length(ebsd.y);
radius = 0.5 * min(max(ebsd.x) - min(ebsd.x), ...
               max(ebsd.y) - min(ebsd.y));
inset = max(numel(ebsd.x) * offset * xres, numel(ebsd.y) * offset * yres);
radius = radius - inset; % inset from the edges of the scan
circ_x = radius * cos(thetas) + xcenter;
circ_y = radius * sin(thetas) + ycenter;
polygon = [circ_x' circ_y'];

[G_N, ~, N_A, N, ~, ~, inside_grains2, edge_grains2, ~] = grainsize_areas_planimetric(ebsd, polygon, varargin{:});

% plotting subfunction
if ismember('PlotResults',varargin)
    plot(ebsd, ebsd.orientations); hold on
    plot(edge_grains2.boundary, 'linewidth', 2, 'lineColor', 'black');
    plot(inside_grains2.boundary, 'linewidth', 3, 'lineColor', 'white');
    plot(polygon(:,1), polygon(:,2), 'k', 'linewidth', 3)

end