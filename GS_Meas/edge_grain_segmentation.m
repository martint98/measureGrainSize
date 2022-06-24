function [inside_grains, edge_grains] = edge_grain_segmentation(ebsd, polygon, varargin)
% Segments inside grains from grains touching the edge of a polygon

% select only ebsd data within polygon
ind = inpolygon(ebsd, polygon);
ebsd = ebsd(ind); 

% detect grains -- note that the variable inside_grains does NOT contain
% only inside grains until the end of the function
[inside_grains,ebsd.grainId] = calcGrains(ebsd('indexed'), 'angle', 5*degree, 'unitcell');

if ismember('exclude_twins', varargin)
    inside_grains = exclude_twins(inside_grains);
end

% delete grains that touch the polygon
outerBoundary_id = any(inside_grains.boundary.grainId==0, 2);
grain_id = inside_grains.boundary(outerBoundary_id).grainId;
edge_grains = inside_grains(grain_id(:,2));
inside_grains(grain_id(:,2)) = [];

end