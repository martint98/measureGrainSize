################################
# Author: Tyler Martin
# Date: 6/14/2022
################################
import warnings
from scipy.io import loadmat
import h5py
import numpy as np
import pandas as pd


class Dstruct:
    def __init__(self, user_data):
        self.data = user_data
        self.euler_angles = user_data['DataContainers']['SyntheticVolumeDataContainer']['CellData']['EulerAngles']
        self.mu = user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']['FeatureSize Distribution']['Average']
        self.sigma = user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']['FeatureSize Distribution']['Standard Deviation']
        self.tmp = user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']['BinNumber']
        self.numbins = len(self.tmp)
        self.minsize = min(self.tmp)
        self.maxsize = max(self.tmp)
        self.meansize = np.exp(self.mu + (0.5 * (self.sigma[0]**2)))
        self.binstepsize = np.mean(self.tmp[2:]-self.tmp[1:-1])
        self.res = user_data['DataContainers']['SyntheticVolumeDataContainer']['_SIMPL_GEOMETRY']['SPACING']
        self.dims = user_data['DataContainers']['SyntheticVolumeDataContainer']['_SIMPL_GEOMETRY']['DIMENSIONS'][:]
        self.mincutoff = (np.log(self.minsize) - self.mu) / self.sigma
        self.maxcutoff = (np.log(self.maxsize) - self.mu) / self.sigma
        self.ngrains = user_data['DataContainers']['SyntheticVolumeDataContainer']['Grain Data'].attrs['TupleDimensions']


def load_d3d_slice(dstruct, sliceindex, plane_normal):
    # shape = size(data);
    shape = np.shape(dstruct.euler_angles)
    # xdim = shape(2);
    # ydim = shape(3);
    # zdim = shape(4); % UNUSED NOW, NEEDED LATER. SEE TO-DO BELOW.
    xdim = shape[0]
    ydim = shape[1]
    # zdim = shape[2]   # UNUSED NOW, NEEDED LATER. SEE TO-DO BELOW.

    # Determine which slice to be analyzed from a certain direction
    # if plane_normal == 'z'
    #     slice = squeeze(data(:, :, :, sliceindex));
    # elseif plane_normal == 'y'
    #     slice = squeeze(data(:, :, sliceindex, :));
    # elseif plane_normal == 'x'
    #     slice = squeeze(data(:, sliceindex, :, :));
    if plane_normal == 'z':
        slice = np.squeeze(dstruct.euler_angles[:][:][sliceindex][:])
    elif plane_normal == 'y':
        slice = np.squeeze(dstruct.euler_angles[:][sliceindex][:][:])
    elif plane_normal == 'x':
        slice = np.squeeze(dstruct.euler_angles[sliceindex][:][:][:])
    else:
        warnings.warn(f"Input plane_normal: {plane_normal} is not valid. Valid inputs are 'x', 'y', and 'z'")

    # TODO: The reshaping needs to be generalized for non-cube datasets!!!

    # Reshape the slice into an n x 3 array of presumed Euler angles
    # Eul = reshape(slice, 3, [])';          # Don't think translation necessary because data already formatted properly

    # Generate an nx1 array of x-coordinates
    # x = repmat(1:xdim, [1, ydim])';
    x = np.zeros((xdim, ydim))
    for i in range(ydim):
        for j in range(xdim):
            x[i, j] = j + 1

    # Generate an nx1 array of y-coordinates
    # y = reshape(repmat(1:ydim, [xdim, 1]), 1, [])';
    y = np.zeros((ydim, xdim))
    for i in range(xdim):
        for j in range(ydim):
            y[i, j] = j + 1

    # Generate an nx1 array of phase indices
    # phase = ones(xdim * ydim, 1);
    phase = np.ones((xdim*ydim, 1))

    # TODO: Convert below code to Python when MTEX functionality is developed in Python
    # %--- Import the slice into MTEX ---%
    # % crystal symmetry
    # CS = {...
    #     'notIndexed',...
    #     crystalSymmetry('m-3m', [3 3 3], 'mineral', 'genericCubic', 'color', [0 0 1]),...
    #     'notIndexed'};
    #
    # % plotting convention
    # setMTEXpref('xAxisDirection','north');
    # setMTEXpref('zAxisDirection','outOfPlane');
    #
    # q = rotation.byEuler(Eul(:,1), Eul(:,2), Eul(:,3));
    # prop.x = x;
    # prop.y = y;
    # ebsd = EBSD(q, phase, CS, prop);
    #
    # return ebsd

if __name__ == '__main__':
    # fpth = "C:\\Users\\marti\\Downloads\\Series1_Structure1_LowRes_1.dream3d"
    fpth = "C:\\Users\\marti\\Downloads\\Supplemental DREAM3D Structures\\Series1_Structure1_HighRes_1.dream3d"

    # Load dream.3d file
    user_data = h5py.File(f'{fpth}', 'r')
    user_dstruct = Dstruct(user_data)

    # Load a slice of the data as MTEX ebsd
    ebsd_path = "C:\\Users\\marti\\Downloads\\GS_Meas\\myEBSD_high_res_1.mat"
    ebsd = loadmat(ebsd_path)
    slice_index = 2
    plane_normal = 'z'
    load_d3d_slice(user_dstruct, slice_index, plane_normal)
    # ebsd = load_d3d_slice(user_data, slice_index, plane_normal)

    # Hard code MTEX ebsd attributes into Python
    ebsd_x_path = "C:\\Users\\marti\\Downloads\\GS_Meas\\myEBSD_high_res_1_x.mat"
    ebsd_y_path = "C:\\Users\\marti\\Downloads\\GS_Meas\\myEBSD_high_res_1_y.mat"
    ebsd_x = loadmat(ebsd_x_path)
    ebsd_y = loadmat(ebsd_y_path)

    # Make all slices have exactly the same dimensions for comparable grain sizes
    # res_adjust = 200.0 / (max(ebsd.x) - min(ebsd.x));
    # ebsd.x = res_adjust * ebsd.x;
    # ebsd.y = res_adjust * ebsd.y;
    res_adjust = 200.0 / (np.max(ebsd_x['ebsd_x']) - np.min(ebsd_x['ebsd_x']))
    ebsd_x = res_adjust * ebsd_x['ebsd_x']
    ebsd_y = res_adjust * ebsd_y['ebsd_y']

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
