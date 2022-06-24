function [G_N, G_A, N_A_counted, N, Abar, areas, inside_grains, edge_grains, polygon] = grainsize_areas_planimetric(ebsd, polygon, varargin)

% Perform segmentation
[inside_grains, edge_grains] = edge_grain_segmentation(ebsd, polygon, varargin{:});

% Calculate number of grains according to ASTM E-112 section 11.1
N_inside = numel(inside_grains.area);

% There is some bug in the edge grain segmentation which results in
% duplicate grains. Selecting only unique centroids fixes this.
u = unique(edge_grains.centroid, 'rows');
N_intercepted = numel(u);

% calculate N_A for grain counting approaches
polygon_area = polyarea(polygon(:,1), polygon(:,2));
N_A_counted = (N_inside + 0.5 * N_intercepted) / polygon_area;

% Note that this can alternatively be calculated as 1/N_A by the standard.
% These values are not equivalent, but should be similar.
areas = inside_grains.area;
Abar = mean(areas);
N = numel(areas);

% Calculate ASTM grain size
G_N = G_numgrain(N_A_counted);
G_A = G_meanbarA(Abar);

end