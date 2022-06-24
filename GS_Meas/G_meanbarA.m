function [G2] = G_meanbarA(bar_A)
% calculating ASTM grain size as a function of the mean cross-sectional
% area of grains not bordering the edge

% bar_A was previously calculated:
    % ar = area(grains)  --area of each individual grain
    % ngrains = numel(ar)  --number of grains detected
    % bar_A = sum(ar)/ngrains  --mean value of grain area

A = 2.0*log2(254.0)+1.0;
B = 1.0/(log10(2.0));

G2 = A - B*log10(bar_A);


end

