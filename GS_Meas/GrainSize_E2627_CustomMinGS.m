function [G_A, Abar, n, N_A_measured, avg_px_per_grain_before_threshold, areas] = GrainSize_E2627_CustomMinGS(ebsd, min_px_per_grain, varargin)
% Area based grain size measurement according to ASTM E2627

% Works with similar outputs to other grain size functions except the
% output argument excluded_grains which is the edge_grans as well as those
% eliminated by the grain size threshold.

% make a rectangle
% this is essentially the same as the rectangle polygon from the Saltikov
% planimetric procedure function, but with null offset. In the Saltikov
% procedure, we want to be able to see that we are drawing a rectangle on
% the figure, so we need to inset it from the edges a bit. For the area
% based method, we want to reject the grains that touch the edges of the
% EBSD scan.
offset = 0.0; % zeroed out
xres = (max(ebsd.x) - min(ebsd.x)) / length(ebsd.x);
yres = (max(ebsd.y) - min(ebsd.y)) / length(ebsd.y);
xinset = numel(ebsd.x) * offset * xres;
yinset = numel(ebsd.y) * offset * yres;

polygon = [min(ebsd.x)+xinset, min(ebsd.y)+yinset;
           min(ebsd.x)+xinset, max(ebsd.y)-yinset;
           max(ebsd.x)-xinset, max(ebsd.y)-yinset;
           max(ebsd.x)-xinset, min(ebsd.y)+yinset];

[~, ~, ~, ~, ~, ~, inside_grains, ~, ~] = grainsize_areas_planimetric(ebsd, polygon, varargin{:});

% ASTM E-2627 takes the mean area of grains with over 100 px each, and
% requires the average grain prior to thresholding has at least 500 px.
% Here we allow an arbitrary number of pixel threshold.
px_area = polyarea(ebsd.unitCell(:,1), ebsd.unitCell(:,2));
threshold = min_px_per_grain * px_area;

% Number of pixels per grain before threshold
avg_px_per_grain_before_threshold = mean(inside_grains.area / px_area);

% Remove grains with fewer pixels than the threshold
excluded_grains = inside_grains(inside_grains.area >= threshold);
inside_grains(inside_grains.area < threshold) = [];
areas = inside_grains.area;
n = length(areas);
Abar = mean(areas);

% Calculate the number of analyzed grains per unit area
% NOTE THAT THIS REFLECTS THE NUMBER OF GRAINS ELIMINATED BY THE THRESHOLD
% This value is potentially useful for assessing differences between
% E112 planimetric measurements and the E2627 standard
analyzed_area = polyarea(polygon(:,1), polygon(:,2)) - sum(excluded_grains.area);
N_A_measured = numel(inside_grains.area) / analyzed_area;

G_A = G_meanbarA(Abar);

% plotting subfunction
if ismember('PlotResults',varargin)
    plot(inside_grains.boundary); hold on
    plot(excluded_grains, 'FaceColor', 'k', 'FaceAlpha', 1.0); hold on
    plot(inside_grains)
end

end