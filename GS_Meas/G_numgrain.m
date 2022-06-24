function G = G_numgrain(N_A)
% ASTM grain size as a function of the number of grains per unit area
% See Planimetric or Jeffries' procedures in ASTM E 112

A = 2.0 * log2(254.0)+1.0;
B = 1.0 / log10(2.0);
G = A + B * log10(N_A);

end

