function [G_L, G_PL, MLI, MIC, grains, intercept_lengths, gb_intersection_coordinates, line_intersection_results, triplept_intersection_coordinates, nlines, total_line_length] = grainsize_linint_random(ebsd, min_intercepts, varargin)

% detect grains
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'), 'angle', 5*degree, 'unitcell');

if ismember('exclude_twins',varargin)
    grains = exclude_twins(grains);
end

% smooth grains
grains = grains.smooth;

% calculate step size
stepsize = 2*abs(ebsd.unitCell(1,1));

%--- Generate a number of random lines ---%
% Lines are added to the slice until at least 50 intercepts are
% calculated.

% Start with one random line generation and add lines until the output
% gives at least 50 intercepts.
nlines = 1; intercept_total = 0;
while intercept_total < min_intercepts
    % Update calculations of intercepts
    %[xyints, xync, linints, length_tot, dist, xycoord] = randlin(n, grains, stepsize);
    [P_L, total_line_length, intercept_lengths, gb_intersection_coordinates, line_intersection_results, triplept_intersection_coordinates] = randlin(ebsd, nlines, grains, stepsize, varargin);
    % Allocate intercept counts
    intercept_count = line_intersection_results(:,5);
    % Total number of intercepts
    intercept_total = sum(intercept_count);
    % Add another random line
    nlines = nlines + 1;
end % end of line addition loop
nlines  = nlines -1;

MLI = mean(intercept_lengths); % Mean lineal intercept
MIC = total_line_length / P_L; % Mean intersection count

G_PL = G_meanintl(MIC);
G_L  = G_meanintl(MLI);

end