clc; fclose('all')

% Startup mtex and return to current directory
mtex_pth = 'C:\Users\marti\Downloads\mtex-5.8.1\mtex-5.8.1';
cd(mtex_pth)
startup_mtex
curdir = 'C:\PyRepo\measureGrainSize\GS_Meas';
cd(curdir)

% Load dream.3d file
% fpth = 'C:\Users\marti\Downloads\Series1_Structure1_LowRes_1.dream3d';
fpth = 'C:\Users\marti\Downloads\Supplemental DREAM3D Structures\Series1_Structure1_HighRes_1.dream3d'
[data, dstruct] = load_d3d_file(fpth);

% Load a slice of the data as MTEX ebsd
sliceindex = 2;
plane_normal = 'z';
ebsd = load_d3d_slice(data, sliceindex, plane_normal);

% Make all slices have exactly the same dimensions for
% comparable grain sizes
res_adjust = 200.0 / (max(ebsd.x) - min(ebsd.x));
ebsd.x = res_adjust * ebsd.x;
ebsd.y = res_adjust * ebsd.y;

% Do some grain size measurements!
%[G_S, N_A_S, n_S] = GrainSize_E112_SaltikovPlanimetric(ebsd);
%[G_J, N_A_J, n_J] = GrainSize_E112_JeffriesPlanimetric(ebsd);
%[G_A1, Abar_A1, n_A1, N_A_measured_A1, avg_px_per_grain_after_threshold, areas_A1] = GrainSize_E2627_AsWritten(ebsd);
%[G_A2, Abar_A2, n_A2, N_A_measured_A2, avg_px_per_grain_before_threshold, areas_A2] = GrainSize_E2627_CustomMinGS(ebsd, 0.0);
%[G_L, lbar, n_L_intercepts, intercept_lengths_L] = GrainSize_E112_HeynRandomLineMLI(ebsd);
%[G_PL, P_L, PL_intersection_count, nlines, Heyn_total_line_length] = GrainSize_E112_HeynRandomLinePL(ebsd);
%[G_Hilliard, hilliardIntCount, hilliard_lbar, hilliardCircumference] = GrainSize_E112_Hilliard(ebsd);
%[G_Abrams, abramsIntCount, abrams_lbar, abramsCircumference] = GrainSize_E112_Abrams(ebsd);
%[G_largestGrain, volFraction] = GrainSize_E930_ALA(ebsd, G_S);
