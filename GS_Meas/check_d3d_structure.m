function [avg_px_per_grain, ngrains] = check_d3d_structure(fpth, randsamples)

ttt = tic;

[data, ~] = load_d3d_file(fpth);

ngrains = zeros(randsamples, 1);
avg_px_per_grain = zeros(randsamples, 1);
for i = 1:randsamples

    rdim = randi(3);
    rslice = randi(size(data, rdim));

    if rdim == 3
        dim = 'z';
    elseif rdim == 2
        dim = 'y';
    elseif rdim ==1
        dim = 'x';
    end

    ebsd = load_d3d_slice(data, rslice, dim);

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

    [~, ~, ~, n, ~, ~, inside_grains, ~, ~] = grainsize_areas_planimetric(ebsd, polygon);

    % ASTM E-2627 takes the mean area of grains with over 100 px each, and
    % requires the average grain prior to thresholding has at least 500 px.
    % Here we allow an arbitrary number of pixel threshold.
    px_area = polyarea(ebsd.unitCell(:,1), ebsd.unitCell(:,2));

    % Number of pixels per grain before threshold
    avg_px_per_grain(i) = mean(inside_grains.area / px_area);
    ngrains(i) = n;

end

avg_px_per_grain = mean(avg_px_per_grain);
ngrains = mean(ngrains);
t = toc(ttt);
fprintf('%6.1f pixels per grain for %i grains from %i random samples in %4.1f seconds.\n', avg_px_per_grain, ngrains, randsamples, t)

end




