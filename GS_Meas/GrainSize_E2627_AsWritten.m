function [G_A, Abar, n, N_A_measured, avg_px_per_grain_before_threshold, areas] = GrainSize_E2627_AsWritten(ebsd, varargin)
% Perform ASTM E2627 measurement as written (minimum grain size 100 px)

min_px_per_grain = 100;
[G_A, Abar, n, N_A_measured, avg_px_per_grain_before_threshold, areas] = GrainSize_E2627_CustomMinGS(ebsd, min_px_per_grain, varargin{:});

end