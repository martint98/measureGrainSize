function [P_L, total_line_length, intercept_lengths, gb_intersection_coordinates, line_intersection_results, triplept_intersection_coordinates] = randlin(ebsd, n, grains, stepsize, varargin)
% Generate random lines on EBSD map and measure intersections and intercept
% lengths for lineal intercept grain size measurements
%
% Input parameters
% ----------------
% n : integer
%                       number of random lines to generate
%
% grains : mtex grain2d object
%
% stepsize : ebsd stepsize
%
% Output parameters
% -----------------
% P_L : scalar
%                       proper intercept count, taking into account ends of
%                       lines and triple points
%
% total_line_length: scalar
%                       sum length of all random lines
%
% intercept_lengths: n x 1 array
%                       lengths between intersections
%
% gb_intersection_coordinates: n x 2 array
%                       column 1: x coordinate of intersection
%                       column 2: y coordinate of intersection
%                       column 3: line number (corresponding to rows in
%                                    line_intersection_results)
%
% line_intersection_results: n x 6 array
%                       column 1: x coordinate start of line
%                       column 2: x coordinate end of line
%                       column 3: y coordinate start of line
%                       column 4: y coordinate end of line
%                       column 5: number of intersections recorded (does
%                                 not take into account whether 
%                                 intersection is at triple point)
%                       column 6: length of line
%
% triplept_intersection_coordinates: n x 2 array
%                       column 1: x coordinate of intersection coincident
%                                 with triple point
%                       column 2: y coordinate of intersection coincident
%                                 with triple point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Have mtex get the coordinates of triple points
tP = grains.triplePoints;
xii = tP.x;
yii = tP.y;
tpoint = [xii, yii];
ntpoints = size(tpoint);
ntpoints = ntpoints(1);

% Get scan dimensions
xdim_max = ceil(max(grains.x)); % maximum x-dimension
ydim_max = ceil(max(grains.y)); % maximum y-dimension
xdim_min = floor(min(grains.x)); % minimum x-dimension
ydim_min = floor(min(grains.y)); % minimum y-dimension

% Begin loop over desired number of random lines
line_intersection_results = []; 
gb_intersection_coordinates = []; % x-y intercepts
total_line_length = 0;
for k = 1:n
    
    % boundaries of the structure
    xdim = [xdim_min, xdim_max]; % x-coordinates of the boundary
    ydim = [ydim_min, ydim_max]; % y-coordinates of the boundary
    
    y = (ydim(2) - ydim(1)) .* rand(2,1) + ydim(1); 
    x = (xdim(2) - xdim(1)) .* rand(2,1) + xdim(1);
    x2 = x(2);
    x1 = x(1);
    y2 = y(2);
    y1 = y(1);
    
    % get slope and intercept of line
    m = (y2 - y1) / (x2 - x1);
    b = y2 - m * x2; % intercept from y=mx+b
    
    % get intersections with bounding box
    yya = m * xdim(1) + b; % value of y at left edge on line
    if yya > ydim(2) % then the x1 coordinate is along top edge
        bbx1 = (ydim(2) - b) / m;
    elseif yya < ydim(1) % then x1 coordinate is along bottom edge
        bbx1 = (ydim(1) - b) / m;
    else % then x coordinate is the left edge
        bbx1 = xdim(1);
    end

    yyb = m * xdim(2) + b; % value of y at right edge on line
    if yyb > ydim(2) % then the x2 coordinate is along top edge
        bbx2 = (ydim(2) - b) / m;
    elseif yyb < ydim(1) % then x2 coordinate is along bottom edge
        bbx2 = (ydim(1) - b) / m;
    else % then x coordinate is the right edge
        bbx2 = xdim(2);
    end

    xxa = (ydim(1) - b) / m; % value of x at y1 on line
    if xxa > xdim(2) % then the y1 coordinate is along right edge
        bby1 = xdim(2) * m + b;
    elseif xxa < xdim(1) % then the y2 coordinate is along left edge
        bby1 = xdim(1) * m + b;
    else % then y coordinate is the bottom edge
        bby1 = ydim(1);
    end
        
    xxb = (ydim(2) - b) / m; % value of x on line at upper edge of bounding box
    if xxb > xdim(2) % then the y2 coordinate is along right edge
        bby2 = xdim(2) * m + b;
    elseif xxb < xdim(1) % then the y2 coordinate is along left edge
        bby2 = xdim(1) * m + b;
    else % it must be the top edge
        bby2 = ydim(2); 
    end
 
    % Collect our line starting and ending points and correct for slope
    offset = 1.0;
    if m>0
        xy1 = [bbx1, bby1] + offset*stepsize;
        xy2 = [bbx2, bby2] - offset*stepsize;
    else
        xy1 = [bbx1 + offset, bby2 - offset*stepsize];
        xy2 = [bbx2 - offset, bby1 + offset*stepsize];
    end
    
    % Have mtex get the intersections
    [xi,yi] = grains.boundary.intersect(xy1,xy2);
   
    % find the number of boundary intersection points
    int_count = sum(~isnan(xi));
    
    % get the x- and y-coordinates of the interceptions
    xx1 = xi(~isnan(xi));
    yy1 = yi(~isnan(yi));
    line_no = k * ones(size(xx1));
    gb_intersection_coordinates = cat(1, gb_intersection_coordinates, [xx1', yy1', line_no']);
    
    % total length of the line
    tot = sqrt((xy2(2) - xy1(2)).^2 + (xy2(1) - xy1(1)).^2);
    total_line_length = total_line_length + tot;
    
    % collate info from individual lines
    line_intersection_results = cat(1, line_intersection_results, [xy1(1), xy1(2), xy2(1), xy2(2), int_count, tot]);
       
end % end of loop over number of lines

% calculate the distance between intersection points and triple points
triplept_intersection_coordinates = [];
tp_thresh = 1.0; % multiples of step size
for m = 1:ntpoints
    % distance in microns:
    dist = sqrt((tpoint(m,1) - gb_intersection_coordinates(1:end,1)).^2 + ...
                (tpoint(m,2) - gb_intersection_coordinates(1:end,2)).^2) * ...
                 tp_thresh * stepsize; 

    %find the distance under threshold and use that as an index into xyints:
    coord = gb_intersection_coordinates(dist<stepsize, :); 
    xcoord = coord(:, 1);
    ycoord = coord(:, 2);
    triplept_intersection_coordinates = cat(1, triplept_intersection_coordinates, [xcoord, ycoord]);
end

% get the count of intersections through the triple points (from xcoord
% and ycoord)
tpcount = numel(xcoord);

% Count the intersections: the ends count as half, hence the -1; 
% add 0.5 counts for each time the line goes through a triple point.
P_L = sum(line_intersection_results(:, 5)) + 0.5 * tpcount - 1; 

% Calculate the intercept lengths
intercept_lengths = sqrt((gb_intersection_coordinates(1:end-1, 1) - ...
                          gb_intersection_coordinates(2:end, 1)).^2 + ...
                         (gb_intersection_coordinates(1:end-1, 2) - ...
                          gb_intersection_coordinates(2:end, 2)).^2);

% plotting subfunction
    if ismember('PlotResults',varargin{:})
        %--- plotting the structure
        plot(ebsd, ebsd.orientations); hold on
        plot(grains.boundary,'LineWidth',1); hold on
        %--- plot triple points
        plot(tP,'color','b','linewidth',2); hold on
        %--- plotting the lines
        plot(ebsd('genericCubic'), ebsd('genericCubic').orientations); hold on
        plot(grains.boundary,'LineWidth',2); hold on
        x = gb_intersection_coordinates(:,1);
        y = gb_intersection_coordinates(:,2); hold on
        for i=1:size(line_intersection_results,1)
            line([line_intersection_results(i,1);line_intersection_results(i,3)], ...
                [line_intersection_results(i,2);line_intersection_results(i,4)], ...
                'linestyle','-','linewidth',4,'color','black')
        end
        %--- plotting the intersections
        hold on
        scatter(x,y,'w','linewidth',2); hold on
        %--- plotting the coordinates considered intersecting a triple
        %point
        xc = triplept_intersection_coordinates(:,1);
        yc = triplept_intersection_coordinates(:,2);
        scatter(xc,yc,'r','linewidth',2)
    end
                      

end % end of the function