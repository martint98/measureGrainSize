function [G_PL, abramsIntCount, abrams_lbar, abramsCircumference_tot] = GrainSize_E112_Abrams(ebsd)
% Abrams Three-Circle Procedure: Three concentric and equally spaced
% circles. 

% Placement of the three-circle test grid should yield 40-100 intersection
% counts per slice tested.

% If the circle intersects a triple point, the count is 2.

% The ratio of circumference is 3:2:1

%--- Populate the grains
[grains, ebsd.grainId] = calcGrains(ebsd('indexed'), 'angle', 2*degree, 'unitcell');
% plot(ebsd, ebsd.orientations); % plot the grains
% hold on
% plot(grains.boundary,'LineWidth',1); % plot the grain boundaries
% hold on 
stepsize = 2*abs(ebsd.unitCell(1,1)); % calculate step size for triple point calculations

% extract triple points
tP = grains.triplePoints;
x_tP = tP.x;
y_tP = tP.y;
tpoint = [x_tP, y_tP];
ntpoints = size(tpoint);
ntpoints = ntpoints(1);

% plot triple points
% plot(tP,'color','b','linewidth',2); hold on


%--- Plotting the largest circle:
        % approximate a circle as a polygon
        offset = 0.02; % 2pct inset from edges
        xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x));
        ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y));
        thetas = 0.0:pi/100.0:2.0*pi;
        xres = 2.0 * xcenter / length(ebsd.x);
        yres = 2.0 * ycenter / length(ebsd.y);
        radius_lg = 0.5 * min(max(ebsd.x) - min(ebsd.x), ...
                       max(ebsd.y) - min(ebsd.y));
        inset = max(numel(ebsd.x) * offset * xres, numel(ebsd.y) * offset * yres);
        radius_lg = radius_lg - inset; % inset from the edges of the scan
        circumference_lg = 2 * pi * radius_lg;
        circ_x_lg = radius_lg * cos(thetas) + xcenter;
        circ_y_lg = radius_lg * sin(thetas) + ycenter;
        polygon_lg = [circ_x_lg' circ_y_lg']; % x/y coords of each line segment

        % plot the largest circle
%         plot(circ_x_lg,circ_y_lg, 'k', 'linewidth', 3)
%         hold on
        
        %--- extract the grain boundary/circle intersection data
        g = size(polygon_lg);
        gg = g(1);

        abrams_intersections_lg = [];
        for n = 1:gg-1
            x_start = polygon_lg(n,1);
            x_end = polygon_lg(n+1,1);
            y_start = polygon_lg(n,2);
            y_end = polygon_lg(n+1,2);
            xy1 = [x_start, y_start];
            xy2 = [x_end, y_end];
            [xi, yi] = grains.boundary.intersect(xy1, xy2);
            x1 = xi(~isnan(xi));
            y1 = yi(~isnan(yi));
            intersect_coords_lg = [x1',y1'];
            %scatter(intersect_coords_lg(:,1),intersect_coords_lg(:,2),'w','linewidth',2)
            num_int_lg(n) = numel(x1);
            abrams_intersections_lg = cat(1, abrams_intersections_lg, intersect_coords_lg);
        end % line intersection loop
        abramsIntCount_lg = sum(num_int_lg)-1;


%--- plotting the medium circle:
        % approximate a circle as a polygon
        % offset = 0.02; % 2pct inset from edges
        xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x));
        ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y));
        thetas = 0.0:pi/100.0:2.0*pi;
        % xres = 2.0 * xcenter / length(ebsd.x);
        % yres = 2.0 * ycenter / length(ebsd.y);
        circumference_med = circumference_lg / 1.5;
        radius_med = circumference_med / (2 * pi);
        % inset = max(numel(ebsd.x) * offset * xres, numel(ebsd.y) * offset * yres);
        radius_med = radius_med - inset; % inset from the edges of the scan
        circ_x_med = radius_med * cos(thetas) + xcenter;
        circ_y_med = radius_med * sin(thetas) + ycenter;
        polygon_med = [circ_x_med' circ_y_med'];

%         plot(circ_x_med,circ_y_med, 'k', 'linewidth', 3)
%         hold on

        g = size(polygon_med);
        gg = g(1);        
        
        abrams_intersections_med = [];
        for n = 1:gg-1
            x_start = polygon_med(n,1);
            x_end = polygon_med(n+1,1);
            y_start = polygon_med(n,2);
            y_end = polygon_med(n+1,2);
            xy1 = [x_start, y_start];
            xy2 = [x_end, y_end];
            [xi, yi] = grains.boundary.intersect(xy1, xy2);
            x1 = xi(~isnan(xi));
            y1 = yi(~isnan(yi));
            intersect_coords_med = [x1',y1'];
            %scatter(intersect_coords_med(:,1),intersect_coords_med(:,2),'w','linewidth',2)
            num_int_med(n) = numel(x1);
            abrams_intersections_med = cat(1, abrams_intersections_med, intersect_coords_med);
        end % line intersection loop
        abramsIntCount_med = sum(num_int_med)-1;

%--- plotting the smallest circle
        % approximate a circle as a polygon
        % offset = 0.02; % 2pct inset from edges
        xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x));
        ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y));
        thetas = 0.0:pi/100.0:2.0*pi;
        % xres = 2.0 * xcenter / length(ebsd.x);
        % yres = 2.0 * ycenter / length(ebsd.y);
        circumference_sm = circumference_lg / 3;
        radius_sm = circumference_sm / (2 * pi);
        % inset = max(numel(ebsd.x) * offset * xres, numel(ebsd.y) * offset * yres);
        radius_sm = radius_sm - inset; % inset from the edges of the scan
        circ_x_sm = radius_sm * cos(thetas) + xcenter;
        circ_y_sm = radius_sm * sin(thetas) + ycenter;
        polygon_sm = [circ_x_sm' circ_y_sm'];

%         plot(circ_x_sm,circ_y_sm, 'k', 'linewidth', 3)
%         hold on
        
        g = size(polygon_sm);
        gg = g(1);

        abrams_intersections_sm = [];
        for n = 1:gg-1
            x_start = polygon_sm(n,1);
            x_end = polygon_sm(n+1,1);
            y_start = polygon_sm(n,2);
            y_end = polygon_sm(n+1,2);
            xy1 = [x_start, y_start];
            xy2 = [x_end, y_end];
            [xi, yi] = grains.boundary.intersect(xy1, xy2);
            x1 = xi(~isnan(xi));
            y1 = yi(~isnan(yi));
            intersect_coords_sm = [x1',y1'];
%             scatter(intersect_coords_sm(:,1),intersect_coords_sm(:,2),'w','linewidth',2)
            num_int_sm(n) = numel(x1);
            abrams_intersections_sm = cat(1, abrams_intersections_sm, intersect_coords_sm);
        end % line intersection loop
        abramsIntCount_sm = sum(num_int_sm)-1;

% Concatenating polygon data
% abramsPolygon = cat(1, abramsPolygon, polygon_sm, polygon_med, polygon_lg);
% g = size(abramsPolygon);
% gg = g(1);

abrams_intersections = cat(1, abrams_intersections_sm, abrams_intersections_med, ...
    abrams_intersections_lg);



% calculate the distance between intersection points and triple points
    triplept_intersection_coordinates = [];
    tp_thresh = 1.0; % multiples of step size
    for m = 1:ntpoints
        % distance in microns:
        dist = sqrt((tpoint(m,1) - abrams_intersections(1:end,1)).^2 + ...
                    (tpoint(m,2) - abrams_intersections(1:end,2)).^2) * ...
                     tp_thresh * stepsize;

        %find the distance under threshold and use that as an index into xyints:
        coord = abrams_intersections(dist<stepsize, :);
        xcoord = coord(:, 1);
        ycoord = coord(:, 2);
        triplept_intersection_coordinates = cat(1, triplept_intersection_coordinates, [xcoord, ycoord]);
    end % triple point distance loop
% get the count of intersections through the triple points (from xcoord
% and ycoord)

xc = triplept_intersection_coordinates(:,1);
yc = triplept_intersection_coordinates(:,2);
% hold on
% scatter(xc,yc,'r','linewidth',2)

abramsTPcount = numel(xc)-1;

abramsIntCount = abramsIntCount_sm + abramsIntCount_med + abramsIntCount_lg;

abramsIntCount = abramsIntCount + abramsTPcount;

% Total line length = total circumference of circles
abramsCircumference_tot = circumference_lg + circumference_med + circumference_sm;

N_L = abramsIntCount / abramsCircumference_tot;
abrams_lbar = 1/N_L;

G_PL = G_meanintl(N_L);

end

