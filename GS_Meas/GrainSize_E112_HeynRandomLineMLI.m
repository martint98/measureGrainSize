function [G_L, lbar, n, intercept_lengths] = GrainSize_E112_HeynRandomLineMLI(ebsd, varargin)

[G_L, ~, ~, ~, ~, intercept_lengths, ~, ~, ~, ~, ~] = grainsize_linint_random(ebsd, 50, varargin{:});

lbar = mean(intercept_lengths);
n = length(intercept_lengths);

end