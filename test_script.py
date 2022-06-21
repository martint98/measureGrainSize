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


class EBSD:
    def __init__(self, x, y):
        self.x = x['ebsd_x']
        self.y = y['ebsd_y']


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

def grainsize_areas_planimetric(polygon):
    # Perform segmentation
    # [inside_grains, edge_grains] = edge_grain_segmentation(ebsd, polygon, varargin{:});

    # Calculate number of grains according to ASTM E-112 section 11.1
    inside_grains_area = 8
    # N_inside = numel(inside_grains.area);
    inside_grains_area = loadmat("C:\\Users\\marti\\Downloads\\GS_Meas\\myEBSD_high_res_1_inside_grains_area.mat")
    N_inside = len(inside_grains_area['inside_grains_area'])

    # There is some bug in the edge grain segmentation which results in
    # duplicate grains. Selecting only unique centroids fixes this.
    edge_grains_centroid = loadmat("C:\\Users\\marti\\Downloads\\GS_Meas\\myEBSD_high_res_1_edge_grains_centroid.mat")
    # u = unique(edge_grains.centroid, 'rows');
    u = np.unique(edge_grains_centroid['edge_grains_centroid'])
    # N_intercepted = numel(u);
    N_intercepted = len(u)

    # Calculate N_A for grain counting approaches
    # polygon_area = polyarea(polygon(:,1), polygon(:,2));
    polygon_area = np.area(polygon)
    print(polygon_area)
    # N_A_counted = (N_inside + 0.5 * N_intercepted) / polygon_area;
    #
    # % Note that this can alternatively be calculated as 1/N_A by the standard.
    # % These values are not equivalent, but should be similar.
    # areas = inside_grains.area;
    # Abar = mean(areas);
    # N = numel(areas);
    #
    # % Calculate ASTM grain size
    # G_N = G_numgrain(N_A_counted);
    # G_A = G_meanbarA(Abar);
    #
    # end

def GrainSize_E112_SaltikovPlanimetric(ebsd):
    # Saltikov's Planimetric: Count of grains in a test rectangle
    # Grains completely enclosed count as 1, those intercepted by the rectangle
    # count by half, and the corners count as one-quarter.

    # Make a rectangle
    offset = 0.02 # 2 percent inset from edges
    xres = (max(ebsd.x) - min(ebsd.x)) / len(ebsd.x)
    yres = (max(ebsd.y) - min(ebsd.y)) / len(ebsd.y)
    # xinset = numel(ebsd.x) * offset * xres;
    # yinset = numel(ebsd.y) * offset * yres;
    xinset = len(ebsd.x) * offset * xres
    yinset = len(ebsd.y) * offset * yres

    # polygon = [min(ebsd.x)+xinset, min(ebsd.y)+yinset;
    #            min(ebsd.x)+xinset, max(ebsd.y)-yinset;
    #            max(ebsd.x)-xinset, max(ebsd.y)-yinset;
    #            max(ebsd.x)-xinset, min(ebsd.y)+yinset];
    polygon = [[(min(ebsd.x) + xinset)[0], (min(ebsd.y) + yinset)[0]],
               [(min(ebsd.x)+xinset)[0], (max(ebsd.y)-yinset)[0]],
               [(max(ebsd.x)-xinset)[0], (max(ebsd.y)-yinset)[0]],
               [(max(ebsd.x)-xinset)[0], (min(ebsd.y)+yinset)[0]]]

    # [G_N, ~, N_A, N, ~, ~, inside_grains2, edge_grains2, ~] = grainsize_areas_planimetric(ebsd, polygon, varargin{:});
    grainsize_areas_planimetric(polygon)
    # if ismember('PlotResults',varargin)
    #     plot(ebsd, ebsd.orientations); hold on
    #     plot(edge_grains2.boundary, 'linewidth', 2, 'lineColor', 'black');
    #     plot(inside_grains2.boundary, 'linewidth', 3, 'lineColor', 'white');
    #     plot(polygon(:,1), polygon(:,2), 'k', 'linewidth', 3)
    # end % parse varargin

if __name__ == '__main__':
    # fpth = "C:\\Users\\marti\\Downloads\\Series1_Structure1_LowRes_1.dream3d"
    fpth = "C:\\Users\\marti\\Downloads\\Supplemental DREAM3D Structures\\Series1_Structure1_HighRes_1.dream3d"

    # Load dream.3d file
    user_data = h5py.File(f'{fpth}', 'r')
    user_dstruct = Dstruct(user_data)

    # Load a slice of the data as MTEX ebsd
    ebsd_path = "C:\\Users\\marti\\Downloads\\GS_Meas\\myEBSD_high_res_1.mat"
    ebsd = loadmat(ebsd_path)
    # Hard code MTEX ebsd attributes into Python
    ebsd_x_path = "C:\\Users\\marti\\Downloads\\GS_Meas\\myEBSD_high_res_1_x.mat"
    ebsd_y_path = "C:\\Users\\marti\\Downloads\\GS_Meas\\myEBSD_high_res_1_y.mat"
    ebsd_x = loadmat(ebsd_x_path)  # Import raw ebsd.x from MTEX as .mat file
    ebsd_y = loadmat(ebsd_y_path)  # Import raw ebsd.y from MTEX as .mat file
    myEBSD = EBSD(ebsd_x, ebsd_y)
    slice_index = 2
    plane_normal = 'z'
    load_d3d_slice(user_dstruct, slice_index, plane_normal)
    # ebsd = load_d3d_slice(user_data, slice_index, plane_normal)

    # Make all slices have exactly the same dimensions for comparable grain sizes
    # res_adjust = 200.0 / (max(ebsd.x) - min(ebsd.x));
    # ebsd.x = res_adjust * ebsd.x;
    # ebsd.y = res_adjust * ebsd.y;
    res_adjust = 200.0 / (np.max(myEBSD.x) - np.min(myEBSD.x))
    # print(res_adjust)
    myEBSD.x = res_adjust * myEBSD.x
    myEBSD.y = res_adjust * myEBSD.y

    # Do some grain size measurements!
    # [G_S, N_A_S, n_S] = GrainSize_E112_SaltikovPlanimetric(ebsd);
    GrainSize_E112_SaltikovPlanimetric(myEBSD)
    # [G_J, N_A_J, n_J] = GrainSize_E112_JeffriesPlanimetric(ebsd);
    # [G_A1, Abar_A1, n_A1, N_A_measured_A1, avg_px_per_grain_after_threshold, areas_A1] = GrainSize_E2627_AsWritten(ebsd);
    # [G_A2, Abar_A2, n_A2, N_A_measured_A2, avg_px_per_grain_before_threshold, areas_A2] = GrainSize_E2627_CustomMinGS(ebsd, 0.0);
    # [G_L, lbar, n_L_intercepts, intercept_lengths_L] = GrainSize_E112_HeynRandomLineMLI(ebsd);
    # [G_PL, P_L, PL_intersection_count, nlines, Heyn_total_line_length] = GrainSize_E112_HeynRandomLinePL(ebsd);
    # [G_Hilliard, hilliardIntCount, hilliard_lbar, hilliardCircumference] = GrainSize_E112_Hilliard(ebsd);
    # [G_Abrams, abramsIntCount, abrams_lbar, abramsCircumference] = GrainSize_E112_Abrams(ebsd);
    # [G_largestGrain, volFraction] = GrainSize_E930_ALA(ebsd, G_S);
