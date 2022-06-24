function [G_PL, hilliardIntCount, hilliard_lbar, circumference] = GrainSize_E112_Hilliard(ebsd, varargin)
%-- Hilliard Single-Circle Procedure: Count of grains intercepting a circle.
% Diameter should never be smaller than the largest observed grains. Test
% circle should be at least 3x the length of the mean lineal intercept.
% Recommended: test conditions that produce ~35 counts per circle.

%-- This function generates a single circle, counts the grains intersecting the
% circle (E112_HilliardGrainCount.m), calculates the mean lineal intercept,
% and plots the data.
%%
% populate grain data
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle', 5*degree,'unitcell');

if ismember('exclude_twins',varargin)
    grains = exclude_twins(grains);
end

stepsize = 2*abs(ebsd.unitCell(1,1));

% extract triple points
tP = grains.triplePoints;
x_tP = tP.x;
y_tP = tP.y;
tpoint = [x_tP, y_tP];
ntpoints = size(tpoint);
ntpoints = ntpoints(1);

%--- generating the circle
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
hilliardPolygon = [circ_x' circ_y'];

%--- generating the points where the line and the grain boundary intersect
hilliard_intersections = [];
for n = 1:length(thetas)-1
    x_start = hilliardPolygon(n,1);
    x_end = hilliardPolygon(n+1,1);
    y_start = hilliardPolygon(n,2);
    y_end = hilliardPolygon(n+1,2);
    xy1 = [x_start, y_start];
    xy2 = [x_end, y_end];
    [xi, yi] = grains.boundary.intersect(xy1, xy2);
    x1 = xi(~isnan(xi));
    y1 = yi(~isnan(yi));
    intersect_coords = [x1',y1'];
    num_int(n) = numel(x1);
    hilliard_intersections = cat(1, hilliard_intersections, intersect_coords);
end % line intersection loop
hilliardIntCount = sum(num_int)-1;
    
% calculate the distance between intersection points and triple points
triplept_intersection_coordinates = [];
tp_thresh = 1.0; % multiples of step size
for m = 1:ntpoints
    % distance in microns:
    dist = sqrt((tpoint(m,1) - hilliard_intersections(1:end,1)).^2 + ...
                (tpoint(m,2) - hilliard_intersections(1:end,2)).^2) * ...
                 tp_thresh * stepsize;

    %find the distance under threshold and use that as an index into xyints:
    coord = hilliard_intersections(dist<stepsize, :);
    xcoord = coord(:, 1);
    ycoord = coord(:, 2);
    triplept_intersection_coordinates = cat(1, triplept_intersection_coordinates, [xcoord, ycoord]);
end % triple point distance loop
% get the count of intersections through the triple points (from xcoord
% and ycoord)
xc = triplept_intersection_coordinates(:,1);
yc = triplept_intersection_coordinates(:,2);
hilliardTPcount = numel(xc)-1;

% % get the count of intersections through the triple points (from xcoord
% % and ycoord)
% tpcount = numel(xcoord);

% Add 0.5 counts for each time the line goes through a triple point
hilliardIntCount = hilliardIntCount + 0.5*hilliardTPcount;

% mean lineal intercept = circumference of circle / number of grains
% intersecting circle
hilliard_lbar = (2*pi*radius) / hilliardIntCount;

G_PL = G_meanintl(hilliardIntCount);

circumference = 2.0 * pi * radius;

%%
% Plotting Subfunction
    if ismember('PlotResults',varargin)
        %--- plot the grains and grain boundaries
        plot(ebsd, ebsd.orientations); hold on
        plot(grains.boundary,'LineWidth',1); hold on
        %--- plot the triple points
        plot(tP,'color','b','linewidth',2); hold on
        %--- plot the circle
        plot(hilliardPolygon(:,1), hilliardPolygon(:,2), ...
            'k', 'linewidth', 3); hold on
        %--- plotting the intersections
        for n = 1:length(thetas)-1
            x_start = hilliardPolygon(n,1);
            x_end = hilliardPolygon(n+1,1);
            y_start = hilliardPolygon(n,2);
            y_end = hilliardPolygon(n+1,2);
            xy1 = [x_start, y_start];
            xy2 = [x_end, y_end];
            [xi, yi] = grains.boundary.intersect(xy1, xy2);
            x1 = xi(~isnan(xi));
            y1 = yi(~isnan(yi));
            intersect_coords = [x1',y1'];
            scatter(intersect_coords(:,1), ...
                intersect_coords(:,2),'w','linewidth',2); hold on
        end % line intersection loop
        scatter(xc,yc,'r','linewidth',2)
        
    end % varargin parse

end % function

