function [G_largestGrain, volFraction] = GrainSize_E930_ALA(ebsd, G2, varargin)
% ASTM Standard E930 Procedure 6.4 (Referee Procedure for Image Analysis)

% G2 ---G number calculated from G_meanbarA.m function

[grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle', 5*degree,'unitcell');

if ismember('exclude_twins',varargin)
    grains = exclude_twins(grains);
end

areas = grains.area;

%--- Calculating the ASTM G number for the largest grain detected
% extract the largest area detected
maxGrainArea = max(areas);

% calculate the ASTM G number for the area in square microns
A = 2.0 * log2(254.0) + 1.0;
B = 1.0 / (log10(2.0));
G_largestGrain = A - B * log10(maxGrainArea);

%---
% G number that represents the average grain size from average area
G_scanAvg = G2;
G_scanAvg = G_scanAvg - 3.0; % subtract 3 ASTM numbers

% convert G_scanAvg back to area
calcArea = 10.0 ^ ((G_scanAvg - A) / -B);

% sum of the areas larger than 'calcArea'
largestGrains = sum(areas(areas>calcArea,:));

% calculating the volume fraction
volFraction = largestGrains / (sum(areas));

% if largest grains make up more than 5% of the area, output a warning
if volFraction > 0.05
    warning('Volume Fraction greater than 0.05')
end

end

