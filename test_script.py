################################
# Author: Tyler Martin
# Date: 6/14/2022
################################

import h5py
import dream3d
import numpy as np


class Dstruct:
    def __init__(self, user_data):
        self.euler_angles = user_data['DataContainers']['SyntheticVolumeDataContainer']['CellData']['EulerAngles']
        self.mu = user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']['FeatureSize Distribution']['Average']
        self.sigma = user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']['FeatureSize Distribution']['Standard Deviation']
        self.numbins = len(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']['BinNumber'])
        self.minsize = min(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']['BinNumber'])
        self.maxsize = max(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']['BinNumber'])
        self.meansize = np.exp(self.mu + (0.5 * (self.sigma**2)))
        self.binstepsize = np.mean(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']['BinNumber'][2:]-user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']['BinNumber'][1:-2])


def load_d3d_file(filepath):
    # Translated from load_d3d_file.m

    # Import h5 data from structure file
    # Read the dataset
    # data = h5read(fpth, '/DataContainers/SyntheticVolumeDataContainer/CellData/EulerAngles');
    data = h5py.File(f'{filepath}', 'r')

    # Extract synthetic structure generation stats
    # dstruct.mu = h5read(fpth, '/DataContainers/StatsGeneratorDataContainer/CellEnsembleData/Statistics/1/FeatureSize Distribution/Average');
    # dstruct.sigma = h5read(fpth, '/DataContainers/StatsGeneratorDataContainer/CellEnsembleData/Statistics/1/FeatureSize Distribution/Standard Deviation');
    # tmp = h5read(fpth, '/DataContainers/StatsGeneratorDataContainer/CellEnsembleData/Statistics/1/BinNumber');
    # dstruct.numbins = length(tmp);
    # dstruct.minsize = min(tmp);
    # dstruct.maxsize = max(tmp);
    # dstruct.meansize = exp(dstruct.mu + 0.5*dstruct.sigma^2);
    # dstruct.binstepsize = mean(tmp(2:end)-tmp(1:end-1));
    # dstruct.res = h5read(fpth, '/DataContainers/SyntheticVolumeDataContainer/_SIMPL_GEOMETRY/SPACING');
    # dstruct.dims = double(h5read(fpth, '/DataContainers/SyntheticVolumeDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));
    # dstruct.mincutoff = (log(dstruct.minsize) - dstruct.mu) / dstruct.sigma;
    # dstruct.maxcutoff = (log(dstruct.maxsize) - dstruct.mu) / dstruct.sigma;
    # dstruct.ngrains = double(h5readatt(fpth, '/DataContainers/SyntheticVolumeDataContainer/Grain Data', 'TupleDimensions'));
    #
    # end

fpth = "C:\\Users\\marti\\Downloads\\Series1_Structure1_LowRes_1.dream3d"

## Startup mtex and return to current directory
# curdir = pwd;
# mtex_pth = "C:\\Users\\marti\\Downloads\\mtex-5.8.1\\mtex-5.8.1"
# cd(mtex_pth)
# startup_mtex
# cd(curdir)

## Load dream.3d file
#[data, dstruct] = load_d3d_file(fpth);
# data, dstruct = load_d3d_file(fpth)
user_data = h5py.File(f'{fpth}', 'r')
user_dstruct = Dstruct(user_data)
# load_d3d_file(fpth)

# Load a slice of the data as MTEX ebsd
# sliceindex = 2;
# plane_normal = 'z';
# ebsd = load_d3d_slice(data, sliceindex, plane_normal);

## Make all slices have exactly the same dimensions for comparable grain sizes
# res_adjust = 200.0 / (max(ebsd.x) - min(ebsd.x));
# ebsd.x = res_adjust * ebsd.x;
# ebsd.y = res_adjust * ebsd.y;

## Do some grain size measurements!
# [G_S, N_A_S, n_S] = GrainSize_E112_SaltikovPlanimetric(ebsd);
# [G_J, N_A_J, n_J] = GrainSize_E112_JeffriesPlanimetric(ebsd);
# [G_A1, Abar_A1, n_A1, N_A_measured_A1, avg_px_per_grain_after_threshold, areas_A1] = GrainSize_E2627_AsWritten(ebsd);
# [G_A2, Abar_A2, n_A2, N_A_measured_A2, avg_px_per_grain_before_threshold, areas_A2] = GrainSize_E2627_CustomMinGS(ebsd, 0.0);
# [G_L, lbar, n_L_intercepts, intercept_lengths_L] = GrainSize_E112_HeynRandomLineMLI(ebsd);
# [G_PL, P_L, PL_intersection_count, nlines, Heyn_total_line_length] = GrainSize_E112_HeynRandomLinePL(ebsd);
# [G_Hilliard, hilliardIntCount, hilliard_lbar, hilliardCircumference] = GrainSize_E112_Hilliard(ebsd);
# [G_Abrams, abramsIntCount, abrams_lbar, abramsCircumference] = GrainSize_E112_Abrams(ebsd);
# [G_largestGrain, volFraction] = GrainSize_E930_ALA(ebsd, G_S);
