function [xbar, s, CI95, RApct] = GrainSize_FieldStats(x_i)
% X_i is an array of the mean values of N_A, P_L, or lbar

% --------------
% Section 15 in ASTM E112
n = length(x_i);
xbar = mean(x_i);
s = sqrt(sum((x_i - xbar).^2) / (n - 1.0));
t = tstat3(n - 1, 0.975, 'inv');
CI95 = t * s / sqrt(n);
RApct = CI95 * 100.0 / xbar;