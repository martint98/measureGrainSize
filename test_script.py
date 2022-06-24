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
        # self.unitcell_x =
        # self.unitcell_y =


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

def G_numgrain(N_A):
    # ASTM grain size as a function of the number of grains per unit area
    # See Planimetric or Jeffries' procedures in ASTM E 112
    A = 2.0 * np.log2(254.0)+1.0
    B = 1.0 / np.log10(2.0)
    G = A + B * np.log10(N_A)

    return G

def G_meanbarA(bar_A):
    # Calculating ASTM grain size as a function of the mean cross-sectional
    # area of grains not bordering the edge

    # TODO: Verify if below note from MATLAB is true for Python code in existence
    # bar_A was previously calculated:
        # ar = area(grains)  --area of each individual grain
        # ngrains = numel(ar)  --number of grains detected
        # bar_A = sum(ar)/ngrains  --mean value of grain area

    A = 2.0*np.log2(254.0)+1.0
    B = 1.0/(np.log10(2.0))

    G2 = A - B*np.log10(bar_A)

    return G2

def grainsize_areas_planimetric(polygon):
    # Perform segmentation
    # [inside_grains, edge_grains] = edge_grain_segmentation(ebsd, polygon, varargin{:});

    # Calculate number of grains according to ASTM E-112 section 11.1
    inside_grains_area = 8
    # N_inside = numel(inside_grains.area);
    inside_grains_area = loadmat("GS_Meas\\myEBSD_high_res_1_inside_grains_area.mat")
    N_inside = len(inside_grains_area['inside_grains_area'])

    # There is some bug in the edge grain segmentation which results in
    # duplicate grains. Selecting only unique centroids fixes this.
    edge_grains_centroid = loadmat("GS_Meas\\myEBSD_high_res_1_edge_grains_centroid.mat")
    # u = unique(edge_grains.centroid, 'rows');
    u = np.unique(edge_grains_centroid['edge_grains_centroid'])
    # N_intercepted = numel(u);
    N_intercepted = len(u)

    # Calculate N_A for grain counting approaches
    # polygon_area = polyarea(polygon(:,1), polygon(:,2));
    # polygon_x = polygon(:, )
    polygon_area = get_polygon_area_shoelace_formula(polygon[:, 0], polygon[:, 1])
    # print(f"polygon_area = {polygon_area}")
    # N_A_counted = (N_inside + 0.5 * N_intercepted) / polygon_area;
    N_A_counted = (N_inside + 0.5 * N_intercepted) / polygon_area

    # Note that this can alternatively be calculated as 1/N_A by the standard.
    # These values are not equivalent, but should be similar.
    # areas = inside_grains.area;
    areas = inside_grains_area['inside_grains_area']
    # Abar = mean(areas);
    Abar = np.mean(areas)
    # N = numel(areas);
    N = len(areas)

    # Calculate ASTM grain size
    # G_N = G_numgrain(N_A_counted);
    G_N = G_numgrain(N_A_counted)
    # G_A = G_meanbarA(Abar);
    G_A = G_meanbarA(Abar)
    return G_N, N_A_counted, N

def get_polygon_area_shoelace_formula(x, y):
    # TODO: Make get_polygon_area_shoelace_formula capable of the below MATLAB polyarea functionality
    # If x and y are vectors of the same length, then polyarea returns the scalar area of the polygon defined by x and y
    # If x and y are matrices of the same size, then polyarea returns a row vector containing the areas of each polygon
    # defined by the columnwise pairs in x and y.
    # If x and y are multidimensional arrays, then polyarea operates along the first dimension whose length is not equal
    # to 1.

    # Assumes x,y points go around the polygon in one direction
    area = abs(sum(i * j for i, j in zip(x, y[1:])) + x[-1] * y[0]
               - sum(i * j for i, j in zip(x[1:], y)) - x[0] * y[-1]) / 2

    return area

def GrainSize_E112_SaltikovPlanimetric(ebsd):
    # Saltikov's Planimetric: Count of grains in a test rectangle
    # Grains completely enclosed count as 1, those intercepted by the rectangle
    # count by half, and the corners count as one-quarter.

    # Make a rectangle
    offset = 0.02  # 2 percent inset from edges
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
    polygon = np.asarray([[(min(ebsd.x) + xinset)[0], (min(ebsd.y) + yinset)[0]],
               [(min(ebsd.x) + xinset)[0], (max(ebsd.y) - yinset)[0]],
               [(max(ebsd.x) - xinset)[0], (max(ebsd.y) - yinset)[0]],
               [(max(ebsd.x) - xinset)[0], (min(ebsd.y) + yinset)[0]]])

    # print(f"Polygon:\n"
    #       f"[[{(min(ebsd.x) + xinset)[0]}, {(min(ebsd.y) + yinset)[0]}],\n"
    #       f"[{(min(ebsd.x)+xinset)[0]}, {(max(ebsd.y)-yinset)[0]}],\n"
    #       f"[{(max(ebsd.x) - xinset)[0]}, {(max(ebsd.y) - yinset)[0]}],\n"
    #       f"[{(max(ebsd.x) - xinset)[0]}, {(min(ebsd.y) + yinset)[0]}]]")

    # [G_N, ~, N_A, N, ~, ~, inside_grains2, edge_grains2, ~] = grainsize_areas_planimetric(ebsd, polygon, varargin{:});
    G_N, N_A, N = grainsize_areas_planimetric(polygon)
    # TODO: Figure out how below MATLAB plotting logic works and implement into Python
    # if ismember('PlotResults',varargin)
    #     plot(ebsd, ebsd.orientations); hold on
    #     plot(edge_grains2.boundary, 'linewidth', 2, 'lineColor', 'black');
    #     plot(inside_grains2.boundary, 'linewidth', 3, 'lineColor', 'white');
    #     plot(polygon(:,1), polygon(:,2), 'k', 'linewidth', 3)
    # end % parse varargin

    return G_N, N_A, N

def GrainSize_E2627_AsWritten(ebsd):
    # function [G_A, Abar, n, N_A_measured, avg_px_per_grain_before_threshold, areas] = GrainSize_E2627_AsWritten(ebsd, varargin)
    # Perform ASTM E2627 measurement as written (minimum grain size 100 px)

    # min_px_per_grain = 100;
    min_px_per_grain = 100
    # [G_A, Abar, n, N_A_measured, avg_px_per_grain_before_threshold, areas] = GrainSize_E2627_CustomMinGS(ebsd, min_px_per_grain, varargin{:});
    G_A, Abar, n, N_A_measured, avg_px_per_grain_before_threshold, areas = GrainSize_E2627_CustomMinGS(ebsd, min_px_per_grain)

    return var1, var2, var3, var4, var5, var6

def GrainSize_E2627_CustomMinGS(ebsd, min_px_per_grain):
    # function [G_A, Abar, n, N_A_measured, avg_px_per_grain_before_threshold, areas] = GrainSize_E2627_CustomMinGS(ebsd, min_px_per_grain, varargin)
    # Area based grain size measurement according to ASTM E2627

    # Works with similar outputs to other grain size functions except the
    # output argument excluded_grains which is the edge_grans as well as those
    # eliminated by the grain size threshold.

    # make a rectangle
    # this is essentially the same as the rectangle polygon from the Saltikov
    # planimetric procedure function, but with null offset. In the Saltikov
    # procedure, we want to be able to see that we are drawing a rectangle on
    # the figure, so we need to inset it from the edges a bit. For the area
    # based method, we want to reject the grains that touch the edges of the
    # EBSD scan.
    # offset = 0.0; % zeroed out
    offset = 0.0  # zeroed out
    # xres = (max(ebsd.x) - min(ebsd.x)) / length(ebsd.x);
    xres = (max((ebsd.x) - min(ebsd.x))) / len(ebsd.x)
    # yres = (max(ebsd.y) - min(ebsd.y)) / length(ebsd.y);
    yres = (max(ebsd.y) - min(ebsd.y)) / len(ebsd.y)
    # xinset = numel(ebsd.x) * offset * xres;
    xinset = len(ebsd.x) * offset * xres
    # yinset = numel(ebsd.y) * offset * yres;
    yinset = len(ebsd.y) * offset * yres

    # polygon = [min(ebsd.x)+xinset, min(ebsd.y)+yinset;
    #            min(ebsd.x)+xinset, max(ebsd.y)-yinset;
    #            max(ebsd.x)-xinset, max(ebsd.y)-yinset;
    #            max(ebsd.x)-xinset, min(ebsd.y)+yinset];

    polygon = np.asarray([[(min(ebsd.x) + xinset)[0], (min(ebsd.y) + yinset)[0]],
                          [(min(ebsd.x) + xinset)[0], (max(ebsd.y) - yinset)[0]],
                          [(max(ebsd.x) - xinset)[0], (max(ebsd.y) - yinset)[0]],
                          [(max(ebsd.x) - xinset)[0], (min(ebsd.y) + yinset)[0]]])

    # [~, ~, ~, ~, ~, ~, inside_grains, ~, ~] = grainsize_areas_planimetric(ebsd, polygon, varargin{:});
    G_N, N_A, N = grainsize_areas_planimetric(polygon)

    # ASTM E-2627 takes the mean area of grains with over 100 px each, and
    # requires the average grain prior to thresholding has at least 500 px.
    # Here we allow an arbitrary number of pixel threshold.
    # px_area = polyarea(ebsd.unitCell(:,1), ebsd.unitCell(:,2));
    px_area = get_polygon_area_shoelace_formula()
    # threshold = min_px_per_grain * px_area;

    # Number of pixels per grain before threshold
    # avg_px_per_grain_before_threshold = mean(inside_grains.area / px_area);

    # Remove grains with fewer pixels than the threshold
    # excluded_grains = inside_grains(inside_grains.area >= threshold);
    # inside_grains(inside_grains.area < threshold) = [];
    # areas = inside_grains.area;
    # n = length(areas);
    # TODO: use max(npArray.shape) if above areas array is not a len of 1 in one dimension
    # Abar = mean(areas);

    # Calculate the number of analyzed grains per unit area
    # NOTE THAT THIS REFLECTS THE NUMBER OF GRAINS ELIMINATED BY THE THRESHOLD
    # This value is potentially useful for assessing differences between
    # E112 planimetric measurements and the E2627 standard
    # analyzed_area = polyarea(polygon(:,1), polygon(:,2)) - sum(excluded_grains.area);
    # N_A_measured = numel(inside_grains.area) / analyzed_area;

    # G_A = G_meanbarA(Abar);

    # plotting subfunction
    # if ismember('PlotResults',varargin)
    #     plot(inside_grains.boundary); hold on
    #     plot(excluded_grains, 'FaceColor', 'k', 'FaceAlpha', 1.0); hold on
    #     plot(inside_grains)

    # return G_A, Abar, n, N_A_measured, avg_px_per_grain_before_threshold, areas

def GrainSize_E112_JeffriesPlanimetric(ebsd):
    # Jeffries' Planimetric: Count of grains in a test circle
    # Grains completely enclosed count as 1, those intercepted by the circle count by half.

    # approximate a circle as a polygon
    # offset = 0.02; % 2pct inset from edges
    offset = 0.02  # 2 percent inset from edges
    # xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x));
    xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x))
    # ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y));
    ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y))
    # thetas = 0.0:pi/100.0:2.0*pi;
    thetas = np.arange(start=0, step=(np.pi/100), stop=(2 * np.pi))
    # xres = 2.0 * xcenter / length(ebsd.x);
    xres = (2 * xcenter)/len(ebsd.x)
    # yres = 2.0 * ycenter / length(ebsd.y);
    yres = (2 * ycenter)/len(ebsd.y)
    # radius = 0.5 * min(max(ebsd.x) - min(ebsd.x), ...
    #                max(ebsd.y) - min(ebsd.y));
    radius = 0.5 * min((max(ebsd.x) - min(ebsd.x)), (max(ebsd.y) - min(ebsd.y)))
    # inset = max(numel(ebsd.x) * offset * xres, numel(ebsd.y) * offset * yres);
    inset = max((len(ebsd.x) * offset * xres), (len(ebsd.y) * offset * yres))
    # radius = radius - inset; % inset from the edges of the scan
    radius = radius - inset  # inset from the edges of the scan
    # circ_x = radius * cos(thetas) + xcenter;
    circ_x = radius * np.cos(thetas) + xcenter
    # circ_y = radius * sin(thetas) + ycenter;
    circ_y = radius * np.sin(thetas) + ycenter
    # polygon = [circ_x' circ_y'];
    for i in range(len(circ_x)):
        if i == 0:
            polygon = np.array([circ_x[i], circ_y[i]])
        else:
            coords = np.array([circ_x[i], circ_y[i]])
            polygon = np.vstack((polygon, coords))

    # [G_N, ~, N_A, N, ~, ~, inside_grains2, edge_grains2, ~] = grainsize_areas_planimetric(ebsd, polygon, varargin{:});
    G_N, N_A, N = grainsize_areas_planimetric(polygon)

    # Plotting subfunction
    # if ismember('PlotResults',varargin)
    #     plot(ebsd, ebsd.orientations); hold on
    #     plot(edge_grains2.boundary, 'linewidth', 2, 'lineColor', 'black');
    #     plot(inside_grains2.boundary, 'linewidth', 3, 'lineColor', 'white');
    #     plot(polygon(:,1), polygon(:,2), 'k', 'linewidth', 3)

    return G_N, N_A, N

if __name__ == '__main__':
    # fpth = "C:\\Users\\marti\\Downloads\\Series1_Structure1_LowRes_1.dream3d"
    fpth = "C:\\Users\\marti\\Downloads\\Supplemental DREAM3D Structures\\Series1_Structure1_HighRes_1.dream3d"

    # Load dream.3d file
    user_data = h5py.File(f'{fpth}', 'r')
    user_dstruct = Dstruct(user_data)

    # Load a slice of the data as MTEX ebsd
    ebsd_path = "GS_Meas\\myEBSD_high_res_1.mat"
    ebsd = loadmat(ebsd_path)
    # Hard code MTEX ebsd attributes into Python
    ebsd_x_path = "GS_Meas\\myEBSD_high_res_1_x.mat"
    ebsd_y_path = "GS_Meas\\myEBSD_high_res_1_y.mat"
    ebsd_x = loadmat(ebsd_x_path)  # Import raw ebsd.x from MTEX as .mat file
    ebsd_y = loadmat(ebsd_y_path)  # Import raw ebsd.y from MTEX as .mat file
    myEBSD = EBSD(ebsd_x, ebsd_y)
    # print(user_data.keys())
    # print(user_data['DataContainerBundles'].keys())
    # print(user_data['DataContainers'].keys())
    # print(user_data['DataContainers']['StatsGeneratorDataContainer'].keys())
    # print(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData'].keys())
    # print(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['CrystalStructures'])   # END
    # print(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['PhaseName'])           # END
    # print(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['PhaseTypes'])          # END
    # print(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['ShapeTypes'])          # END
    # print(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics'].keys())
    # print(user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1'].keys())
    # print(user_data['DataContainers']['StatsGeneratorDataContainer']['_SIMPL_GEOMETRY'])                         # END
    # print(user_data['DataContainers']['SyntheticVolumeDataContainer'].keys())
    # print(user_data['DataContainers']['SyntheticVolumeDataContainer']['CellData'].keys())
    # print(user_data['DataContainers']['SyntheticVolumeDataContainer']['CellEnsembleData'].keys())
    # print(user_data['DataContainers']['SyntheticVolumeDataContainer']['Grain Data'].keys())
    # print(user_data['DataContainers']['SyntheticVolumeDataContainer']['_SIMPL_GEOMETRY'].keys())
    # print(user_data['Pipeline'].keys())
    # print(user_data['Pipeline']['Pipeline'])    # END
    # print(f"Max x: {max(myEBSD.x)}\nMax y: {max(myEBSD.y)}\nMin x: {min(myEBSD.x)}\nMin y: {min(myEBSD.y)}")
    slice_index = 2
    plane_normal = 'z'
    load_d3d_slice(user_dstruct, slice_index, plane_normal)
    # ebsd = load_d3d_slice(user_data, slice_index, plane_normal)

    # Make all slices have exactly the same dimensions for comparable grain sizes
    # res_adjust = 200.0 / (max(ebsd.x) - min(ebsd.x));
    # ebsd.x = res_adjust * ebsd.x;
    # ebsd.y = res_adjust * ebsd.y;
    res_adjust = 200.0 / (max(myEBSD.x) - min(myEBSD.x))
    # print(res_adjust)
    myEBSD.x = res_adjust * myEBSD.x
    myEBSD.y = res_adjust * myEBSD.y
    # print(f"res_adjust: {res_adjust}")
    # print(f"Adjusted:\nMax x: {max(myEBSD.x)}\nMax y: {max(myEBSD.y)}\nMin x: {min(myEBSD.x)}\nMin y: {min(myEBSD.y)}")

    # Do some grain size measurements!
    # [G_S, N_A_S, n_S] = GrainSize_E112_SaltikovPlanimetric(ebsd);
    G_N, N_A, N = GrainSize_E112_SaltikovPlanimetric(myEBSD)
    # [G_J, N_A_J, n_J] = GrainSize_E112_JeffriesPlanimetric(ebsd);
    G_N, N_A, N = GrainSize_E112_JeffriesPlanimetric(myEBSD)
    # [G_A1, Abar_A1, n_A1, N_A_measured_A1, avg_px_per_grain_after_threshold, areas_A1] = GrainSize_E2627_AsWritten(ebsd);
    # var1, var2, var3, var4, var5, var6 = GrainSize_E2627_AsWritten(myEBSD)
    # [G_A2, Abar_A2, n_A2, N_A_measured_A2, avg_px_per_grain_before_threshold, areas_A2] = GrainSize_E2627_CustomMinGS(ebsd, 0.0);
    # [G_L, lbar, n_L_intercepts, intercept_lengths_L] = GrainSize_E112_HeynRandomLineMLI(ebsd);
    # [G_PL, P_L, PL_intersection_count, nlines, Heyn_total_line_length] = GrainSize_E112_HeynRandomLinePL(ebsd);
    # [G_Hilliard, hilliardIntCount, hilliard_lbar, hilliardCircumference] = GrainSize_E112_Hilliard(ebsd);
    # [G_Abrams, abramsIntCount, abrams_lbar, abramsCircumference] = GrainSize_E112_Abrams(ebsd);
    # [G_largestGrain, volFraction] = GrainSize_E930_ALA(ebsd, G_S);
