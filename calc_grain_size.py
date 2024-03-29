################################
# Author: Tyler Martin
# Date: 6/14/2022
################################
import math
import warnings
from scipy.io import loadmat
import h5py
import numpy as np
from matplotlib import path


class Dstruct:
    def __init__(self, user_data):
        self.data = user_data
        self.euler_angles = user_data['DataContainers']['SyntheticVolumeDataContainer']['CellData']['EulerAngles']
        self.temp = user_data['DataContainers']['StatsGeneratorDataContainer']['CellEnsembleData']['Statistics']['1']
        self.mu = self.temp['FeatureSize Distribution']['Average']
        self.sigma = self.temp['FeatureSize Distribution']['Standard Deviation']
        self.tmp = self.temp['BinNumber']
        self.numbins = len(self.tmp)
        self.minsize = min(self.tmp)
        self.maxsize = max(self.tmp)
        self.meansize = np.exp(self.mu + (0.5 * (self.sigma[0]**2)))
        self.binstepsize = np.mean(self.tmp[2:]-self.tmp[1:-1])
        self.res = user_data['DataContainers']['SyntheticVolumeDataContainer']['_SIMPL_GEOMETRY']['SPACING']
        self.dims = user_data['DataContainers']['SyntheticVolumeDataContainer']['_SIMPL_GEOMETRY']['DIMENSIONS'][:]
        self.mincutoff = (np.log(self.minsize) - self.mu) / self.sigma
        self.maxcutoff = (np.log(self.maxsize) - self.mu) / self.sigma
        self.grain_data = user_data['DataContainers']['SyntheticVolumeDataContainer']['Grain Data']
        self.ngrains = self.grain_data.attrs['TupleDimensions']


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

    # Note from MATLAB (keeping for now)
    # bar_A was previously calculated:
    # ar = area(grains)  --area of each individual grain
    # ngrains = numel(ar)  --number of grains detected
    # bar_A = sum(ar)/ngrains  --mean value of grain area

    A = 2.0*np.log2(254.0)+1.0
    B = 1.0/(np.log10(2.0))
    G2 = A - B*np.log10(bar_A)

    return G2


def G_meanintl(u):
    # Calculate the ASTM grain size as a function of mean intercept length

    # mean intercept length was previously calculated:
    # Calculating the average intercept length
    # ints = xync(:,5)  --extracting number of intercept lines
    # z = sum(ints)  --total number of intercept lines
    # intl = linints(:,1)  --extracting the intercept lengths
    # q = sum(intl)  --total of intercept lengths
    # Calculating the mean intercept length:
    # u = q / z  --total of int lengths / total number of intercepts

    A = 2.0 * np.log2(320.0)
    B = 2.0 / np.log10(2.0)

    G = A - B * np.log10(u)

    return G


def edge_grain_segmentation(ebsd, polygon):
    # Segments inside grains from grains touching the edge of a polygon

    # Select only ebsd data within polygon
    p = path.Path(polygon)
    ind = [i for i in ebsd if p.contains_point(ebsd[i])]        # Indices of EBSD points within polygon of interest
    ebsd = ebsd[ind]

    # Detect grains -- note that the variable inside_grains does NOT contain
    # only inside grains until the end of the function
    degree = 1
    # inside_grains, ebsd.grainId = calcGrains(ebsd('indexed'), 'angle', 5 * degree, 'unitcell')

    # if ismember('exclude_twins'):
    #     inside_grains = exclude_twins(inside_grains)

    # Delete grains that touch the polygon
    # outerBoundary_id = np.any(inside_grains.boundary.grainId==0, 2)
    # grain_id = inside_grains.boundary(outerBoundary_id).grainId
    # edge_grains = inside_grains[grain_id[:, 2]]
    # inside_grains[grain_id[:, 2]] = []

    # return inside_grains, edge_grains


def grainsize_areas_planimetric(polygon_x, polygon_y, measurement_type: str):
    # Perform segmentation
    # [inside_grains, edge_grains] = edge_grain_segmentation(ebsd, polygon, varargin{:});

    # Calculate number of grains according to ASTM E-112 section 11.1
    if measurement_type == 'Saltikov' or 'ALA':
        inside_grains_area = loadmat("GS_Meas\\myEBSD_high_res_1_inside_grains_area.mat")  # Saltikov and ALA
        inside_grains_area = inside_grains_area['inside_grains_area']                      # Saltikov and ALA
    if measurement_type == 'Jeffries':
        inside_grains_area = loadmat("GS_Meas\\myEBSD_high_res_1_area_inside_grains.mat")  # Jeffries
        inside_grains_area = inside_grains_area['areaInsideGrains']                        # Jeffries
    N_inside = len(inside_grains_area)

    # There is some bug in the edge grain segmentation which results in
    # duplicate grains. Selecting only unique centroids fixes this.
    # u = unique(edge_grains.centroid, 'rows');
    if measurement_type == 'Saltikov' or 'ALA':
        edge_grains_centroid = loadmat("GS_Meas\\myEBSD_high_res_1_edge_grains_centroid.mat")         # Saltikov and ALA
        u = np.unique(edge_grains_centroid['edge_grains_centroid'])                                   # Saltikov and ALA
    if measurement_type == 'Jeffries':
        edge_grains_centroid = loadmat("GS_Meas\\myEBSD_high_res_1_jeffries_edge_grains_centroid.mat")  # Jeffriees
        u = np.unique(edge_grains_centroid['jeffries_edge_grain_centroids'])                            # Jeffries
    N_intercepted = len(u)

    # Calculate N_A for grain counting approaches
    polygon_area = get_polygon_area_shoelace_formula(polygon_x, polygon_y)
    N_A_counted = (N_inside + 0.5 * N_intercepted) / polygon_area

    # Note that this can alternatively be calculated as 1/N_A by the standard.
    # These values are not equivalent, but should be similar.
    # areas = inside_grains.area;
    areas = inside_grains_area
    Abar = np.mean(areas)
    N = len(areas)

    # Calculate ASTM grain size
    G_N = G_numgrain(N_A_counted)
    G_A = G_meanbarA(Abar)
    return G_N, N_A_counted, N


def get_polygon_area_shoelace_formula(x, y):
    """Calculate area of polygon
    :param x: numpy array of polygon vertices x values
    :param y: numpy array of polygon vertices y values
    :return: area of the polygon
    """
    # Assumes x,y points go around the polygon in one direction
    area = abs(sum(i * j for i, j in zip(x, y[1:])) + x[-1] * y[0]
               - sum(i * j for i, j in zip(x[1:], y)) - x[0] * y[-1]) / 2
    return area


def grainsize_linint_random(ebsd, min_intercepts: int = 50):
    # detect grains
    # [grains,ebsd.grainId] = calcGrains(ebsd('indexed'), 'angle', 5*degree, 'unitcell');

    # TODO: Implement below functionality as it doesn't appear used at the moment
    # if ismember('exclude_twins',varargin)
    #     grains = exclude_twins(grains);

    # smooth grains
    # grains = grains.smooth;
    # grains = smooth(grains, 5)
    grains = loadmat("GS_Meas\\myEBSD_high_res_1_smoothed_grains.mat")          # HeynRandomLineMLI and PL

    # calculate step size
    # stepsize = 2*abs(ebsd.unitCell(1,1));
    unitCell = loadmat("GS_Meas\\myEBSD_high_res_1_unit_cell.mat")              # HeynRandomLineMLI and PL
    unitCell = unitCell['myUnitCell']                                           # HeynRandomLineMLI and PL
    stepsize = 2 * abs(unitCell[0][0])

    # Generate a number of random lines
    # Lines are added to the slice until at least 50 intercepts are
    # calculated.

    # Start with one random line generation and add lines until the output
    # gives at least 50 intercepts.
    # nlines = 1; intercept_total = 0;
    nlines = 1
    intercept_total = 0

    # while intercept_total < min_intercepts
    # # while intercept_total < min_intercepts:
    #      # Update calculations of intercepts
    #      %[xyints, xync, linints, length_tot, dist, xycoord] = randlin(n, grains, stepsize);
    #      P_L, total_line_length, intercept_lengths, gb_intersection_coordinates, line_intersection_results,
    #       triplept_intersection_coordinates = randlin(ebsd, nlines, grains, stepsize, varargin);
    #      TODO: Finish translation of below (skipping for now due to usage of MTEX within randlin func in this
    #       while loop)
    #      P_L, total_line_length, intercept_lengths, gb_intersection_coordinates, line_intersection_results, \
    #       triplept_intersection_coordinates = randlin(ebsd, nlines, grains, stepsize)
    #      Allocate intercept counts
    #      intercept_count = line_intersection_results(:,5);
    #      Total number of intercepts
    #      intercept_total = sum(intercept_count);
    #      Add another random line
    #      nlines = nlines + 1;

    # Hard code loop results for HeynRandomLineMLI and PL
    P_L = 51
    total_line_length = 1081.373455580337
    intercept_lengths = loadmat("GS_Meas\\linint_random_while_loop_output\\intercept_lengths.mat")
    intercept_lengths = intercept_lengths["intercept_lengths"]
    gb_intersection_coordinates = loadmat("GS_Meas\\linint_random_while_loop_output\\gb_intersection_coordinates.mat")
    gb_intersection_coordinates = gb_intersection_coordinates["gb_intersection_coordinates"]
    line_intersection_results = loadmat("GS_Meas\\linint_random_while_loop_output\\line_intersection_coordinates.mat")
    line_intersection_results = line_intersection_results["line_intersection_results"]
    triplept_intersection_coordinates = loadmat("GS_Meas\\linint_random_while_loop_output\\"
                                                "triplept_intersection_coordinates.mat")
    triplept_intersection_coordinates = triplept_intersection_coordinates["triplept_intersection_coordinates"]
    nlines = 6
    intercept_count = line_intersection_results[:, 4]
    intercept_total = sum(intercept_count)

    nlines = nlines - 1

    MLI = np.mean(intercept_lengths)  # Mean lineal intercept
    MIC = total_line_length / P_L  # Mean intersection count

    G_PL = G_meanintl(MIC)
    G_L = G_meanintl(MLI)

    return G_L, G_PL, MLI, MIC, grains, intercept_lengths, gb_intersection_coordinates, line_intersection_results, \
        triplept_intersection_coordinates, nlines, total_line_length


def smooth(a, WSZ):
    out0 = np.convolve(a, np.ones(WSZ, dtype=int), 'valid')/WSZ
    r = np.arange(1, WSZ-1, 2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]

    return np.concatenate((start, out0, stop))


def randlin(ebsd, n, grains, stepsize):
    # % Generate random lines on EBSD map and measure intersections and intercept
    # % lengths for lineal intercept grain size measurements
    # %
    # % Input parameters
    # % ----------------
    # % n : integer
    # %                       number of random lines to generate
    # %
    # % grains : mtex grain2d object
    # %
    # % stepsize : ebsd stepsize
    # %
    # % Output parameters
    # % -----------------
    # % P_L : scalar
    # %                       proper intercept count, taking into account ends of
    # %                       lines and triple points
    # %
    # % total_line_length: scalar
    # %                       sum length of all random lines
    # %
    # % intercept_lengths: n x 1 array
    # %                       lengths between intersections
    # %
    # % gb_intersection_coordinates: n x 2 array
    # %                       column 1: x coordinate of intersection
    # %                       column 2: y coordinate of intersection
    # %                       column 3: line number (corresponding to rows in
    # %                                    line_intersection_results)
    # %
    # % line_intersection_results: n x 6 array
    # %                       column 1: x coordinate start of line
    # %                       column 2: x coordinate end of line
    # %                       column 3: y coordinate start of line
    # %                       column 4: y coordinate end of line
    # %                       column 5: number of intersections recorded (does
    # %                                 not take into account whether
    # %                                 intersection is at triple point)
    # %                       column 6: length of line
    # %
    # % triplept_intersection_coordinates: n x 2 array
    # %                       column 1: x coordinate of intersection coincident
    # %                                 with triple point
    # %                       column 2: y coordinate of intersection coincident
    # %                                 with triple point
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Have mtex get the coordinates of triple points
    # tP = grains.triplePoints;
    x_tP = loadmat("GS_Meas\\myEBSD_high_res_1_x_tp.mat")
    x_tP = x_tP['x_tP']
    y_tP = loadmat("GS_Meas\\myEBSD_high_res_1_y_tp.mat")
    y_tP = y_tP['y_tP']
    tpoint = np.hstack((x_tP, y_tP))
    ntpoints = np.shape(tpoint)
    ntpoints = ntpoints[0]

    # Get scan dimensions
    # xdim_max = ceil(max(grains.x)); % maximum x-dimension
    xdim_max = math.ceil(max(ebsd.x))  # Maximum x-dimension
    # ydim_max = ceil(max(grains.y)); % maximum y-dimension
    ydim_max = math.ceil(max(ebsd.y))  # Maximum y-dimension
    # xdim_min = floor(min(grains.x)); % minimum x-dimension
    xdim_min = math.floor(min(ebsd.x))  # Minimum y-dimension
    # ydim_min = floor(min(grains.y)); % minimum y-dimension
    ydim_min = math.floor(min(ebsd.y))  # Minimum y-dimension

    # Begin loop over desired number of random lines
    line_intersection_results = []
    gb_intersection_coordinates = []  # x-y intercepts
    total_line_length = 0

    for k in range(n):
        xdim = np.array([xdim_min, xdim_max])  # x-coordinates of the boundary
        ydim = np.array([ydim_min, ydim_max])  # y-coordinates of the boundary
        y = (ydim_max - ydim_min) * np.random.randint(0, (2, 1)) + ydim_min
        x = (xdim_max - xdim_min) * np.random.randint(0, (2, 1)) + xdim_min
        # Get slope and intercept of line
        m = (ydim_max - ydim_min) / (xdim_max - xdim_min)
        b = ydim_max - m + xdim_max  # Intercept from y=mx+b
        # Get intersections with bounding box
        yya = m * xdim_min + b  # Value of y at left edge on line
        if yya > ydim_max:  # Then the xdim_min coordinate is along the top edge
            bbx1 = (ydim_max - b) / m
        elif yya < ydim_min:  # Then xdim_min coordinate is along bottom edge
            bbx1 = (ydim_min - b) / m
        else:  # Then x coordinate is the left edge
            bbx1 = xdim_min
        yyb = m * xdim_max + b  # Value of y at right edge on line
        if yyb > ydim_max:  # Then the xdim_max coordinate is along top edge
            bbx2 = (ydim_max - b) / m
        elif yyb < ydim_min:  # Then xdim_max coordinate is along bottom edge
            bbx2 = (ydim_min - b) / m
        else:  # Then x coordinate is the right edge
            bbx2 = xdim_max
        xxa = (ydim_min - b) / m  # Value of x at y1 on line
        if xxa > xdim_max:  # Then the y1 coordinate is along right edge
            bby1 = xdim_max * m + b
        elif xxa < xdim_min:  # Then the y2 coordinate is along left edge
            bby1 = xdim_min * m + b
        else:  # Then y coordinate is the bottom edge
            bby1 = ydim_min
        xxb = (ydim_max - b) / m  # Value of x on line at upper edge of bounding box
        if xxb > xdim_max:  # Then the y2 coordinate is along left edge
            bby2 = xdim_max * m + b
        elif xxb < xdim_min:  # Then the y2 coordinate is along right edge
            bby2 = xdim_min * m + b
        else:  # It must be top edge
            bby2 = ydim_max
        # Collect our line starting and ending points and correct for slope
        offset = 1.0
        if m > 0:
            xy1 = np.array([bbx1, bby1]) + offset * stepsize
            xy2 = np.array([bbx2, bby2]) - offset * stepsize
        else:
            xy1 = np.array([(bbx1 + offset), (bby2 - offset * stepsize)])
            xy2 = np.array([(bbx2 - offset), (bby1 + offset * stepsize)])
        # Have mtex get the intersections
        # TODO: Create implementation of below
        # [xi,yi] = grains.boundary.intersect(xy1,xy2);
        xi = []
        yi = []
        # Find the number of boundary intersection points
        int_count = sum(np.logical_not(np.isnan(xi)))
        # Get the x- and y-coordinates of the interceptions
        xx1 = xi[np.logical_not(np.isnan(xi))]
        yy1 = yi[np.logical_not(np.isnan(yi))]
        line_no = k * np.ones((1, len(xx1)))
        if len(gb_intersection_coordinates) == 0:
            gb_intersection_coordinates = np.array([xx1.T.conj(), yy1.T.conj(), line_no.T.conj()])
        else:
            gb_intersection_coordinates = np.vstack([gb_intersection_coordinates, np.array([xx1.T.conj(), yy1.T.conj(),
                                                                                            line_no.T.conj()])])
        # Total length of the line
        tot = np.sqrt((xy2[1] - xy1[1]) ** 2 + (xy2[0] - xy1[0]) ** 2)
        total_line_length = total_line_length + tot
        # Collate info from individual lines
        if len(line_intersection_results) == 0:
            line_intersection_results = np.array([xy1[0], xy1[1], xy2[0], xy2[1], int_count, tot])
        else:
            line_intersection_results = np.vstack([line_intersection_results, np.array([xy1[0], xy1[1], xy2[0], xy2[1],
                                                                                        int_count, tot])])

    # Calculate the distance between intersection points and triple points
    triplept_intersection_coordinates = []
    tp_thresh = 1.0  # Multiples of step size
    xcoord = 0
    for m in range(ntpoints):
        # Distance in microns
        dist = np.sqrt(((tpoint[m, 0] - gb_intersection_coordinates[:, 0]) ** 2) +
                       ((tpoint[m, 1] - gb_intersection_coordinates[:, 1]) ** 2)) * tp_thresh * stepsize

        # Find the distance under threshold and use that as an index into xyints:
        current_coord = gb_intersection_coordinates[dist < stepsize]
        if len(current_coord) != 0:
            for i in range(len(current_coord)):
                coord = current_coord[i]
                xcoord = coord[:][0]
                ycoord = coord[:][1]
                gb_intersection_coordinates = np.hstack((xcoord, ycoord))
    # Get the count of intersections through the triple points (from xcoord and ycoord)
    tpcount = len(xcoord)

    # Count the intersections: the ends count as half, hence the -1
    # Add 0.5 counts for each time the line goes through a triple point.
    P_L = sum(line_intersection_results[:, 4]) + 0.5 * tpcount - 1

    # Calculate the intercept lengths
    intercept_lengths = np.sqrt((gb_intersection_coordinates[:-2, 0] - gb_intersection_coordinates[1:, 0]) ** 2 +
                                (gb_intersection_coordinates[:-2, 1] - gb_intersection_coordinates[1:, 1]) ** 2)

    # % plotting subfunction
    #     if ismember('PlotResults',varargin{:})
    #         %--- plotting the structure
    #         plot(ebsd, ebsd.orientations); hold on
    #         plot(grains.boundary,'LineWidth',1); hold on
    #         %--- plot triple points
    #         plot(tP,'color','b','linewidth',2); hold on
    #         %--- plotting the lines
    #         plot(ebsd('genericCubic'), ebsd('genericCubic').orientations); hold on
    #         plot(grains.boundary,'LineWidth',2); hold on
    #         x = gb_intersection_coordinates(:,1);
    #         y = gb_intersection_coordinates(:,2); hold on
    #         for i=1:size(line_intersection_results,1)
    #             line([line_intersection_results(i,1);line_intersection_results(i,3)], ...
    #                 [line_intersection_results(i,2);line_intersection_results(i,4)], ...
    #                 'linestyle','-','linewidth',4,'color','black')
    #         end
    #         %--- plotting the intersections
    #         hold on
    #         scatter(x,y,'w','linewidth',2); hold on
    #         %--- plotting the coordinates considered intersecting a triple
    #         %point
    #         xc = triplept_intersection_coordinates(:,1);
    #         yc = triplept_intersection_coordinates(:,2);
    #         scatter(xc,yc,'r','linewidth',2)
    #     end

    return P_L, total_line_length, intercept_lengths, gb_intersection_coordinates, line_intersection_results, \
        triplept_intersection_coordinates


def GrainSize_E112_SaltikovPlanimetric(ebsd):
    # Saltikov's Planimetric: Count of grains in a test rectangle
    # Grains completely enclosed count as 1, those intercepted by the rectangle
    # count by half, and the corners count as one-quarter.

    # Make a rectangle
    offset = 0.02  # 2 percent inset from edges
    xres = (max(ebsd.x) - min(ebsd.x)) / len(ebsd.x)
    yres = (max(ebsd.y) - min(ebsd.y)) / len(ebsd.y)
    xinset = len(ebsd.x) * offset * xres
    yinset = len(ebsd.y) * offset * yres

    polygon = np.asarray([[(min(ebsd.x) + xinset)[0], (min(ebsd.y) + yinset)[0]],
                          [(min(ebsd.x) + xinset)[0], (max(ebsd.y) - yinset)[0]],
                          [(max(ebsd.x) - xinset)[0], (max(ebsd.y) - yinset)[0]],
                          [(max(ebsd.x) - xinset)[0], (min(ebsd.y) + yinset)[0]]])

    # [G_N, ~, N_A, N, ~, ~, inside_grains2, edge_grains2, ~] = grainsize_areas_planimetric(ebsd, polygon, varargin{:});
    G_N, N_A, N = grainsize_areas_planimetric(polygon[:, 0], polygon[:, 1], 'Saltikov')

    # TODO: Figure out how below MATLAB plotting logic works and implement into Python
    # if ismember('PlotResults',varargin)
    #     plot(ebsd, ebsd.orientations); hold on
    #     plot(edge_grains2.boundary, 'linewidth', 2, 'lineColor', 'black');
    #     plot(inside_grains2.boundary, 'linewidth', 3, 'lineColor', 'white');
    #     plot(polygon(:,1), polygon(:,2), 'k', 'linewidth', 3)
    # end % parse varargin

    return G_N, N_A, N


def GrainSize_E112_JeffriesPlanimetric(ebsd):
    # Jeffries' Planimetric: Count of grains in a test circle
    # Grains completely enclosed count as 1, those intercepted by the circle count by half.

    # approximate a circle as a polygon
    offset = 0.02  # 2 percent inset from edges
    xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x))
    ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y))
    thetas = np.arange(start=0, step=(np.pi/100), stop=((2 * np.pi) + (np.pi/100)))
    xres = (2 * xcenter)/len(ebsd.x)
    yres = (2 * ycenter)/len(ebsd.y)
    radius = 0.5 * min((max(ebsd.x) - min(ebsd.x)), (max(ebsd.y) - min(ebsd.y)))
    inset = max((len(ebsd.x) * offset * xres), (len(ebsd.y) * offset * yres))
    radius = radius - inset  # inset from the edges of the scan
    circ_x = radius * np.cos(thetas) + xcenter
    circ_y = radius * np.sin(thetas) + ycenter

    # [G_N, ~, N_A, N, ~, ~, inside_grains2, edge_grains2, ~] = grainsize_areas_planimetric(ebsd, polygon, varargin{:});
    G_N, N_A, N = grainsize_areas_planimetric(circ_x, circ_y, 'Jeffries')

    # Plotting subfunction
    # if ismember('PlotResults',varargin)
    #     plot(ebsd, ebsd.orientations); hold on
    #     plot(edge_grains2.boundary, 'linewidth', 2, 'lineColor', 'black');
    #     plot(inside_grains2.boundary, 'linewidth', 3, 'lineColor', 'white');
    #     plot(polygon(:,1), polygon(:,2), 'k', 'linewidth', 3)

    return G_N, N_A, N


def GrainSize_E2627_AsWritten(ebsd):
    # Perform ASTM E2627 measurement as written (minimum grain size 100 px)

    min_px_per_grain = 100
    G_A, Abar, n, N_A_measured, avg_px_per_grain_before_threshold, areas = \
        GrainSize_E2627_CustomMinGS(ebsd, min_px_per_grain, 'AsWritten')

    return G_A, Abar, n, N_A_measured, avg_px_per_grain_before_threshold, areas


def GrainSize_E2627_CustomMinGS(ebsd, min_px_per_grain: int = 0, measurement_type: str = 'CustomMinGS'):
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
    offset = 0.0  # zeroed out
    xres = max(ebsd.x) - min(ebsd.x) / len(ebsd.x)
    yres = max(ebsd.y) - min(ebsd.y) / len(ebsd.y)
    xinset = len(ebsd.x) * offset * xres
    yinset = len(ebsd.y) * offset * yres

    polygon = np.asarray([[(min(ebsd.x) + xinset)[0], (min(ebsd.y) + yinset)[0]],
                          [(min(ebsd.x) + xinset)[0], (max(ebsd.y) - yinset)[0]],
                          [(max(ebsd.x) - xinset)[0], (max(ebsd.y) - yinset)[0]],
                          [(max(ebsd.x) - xinset)[0], (min(ebsd.y) + yinset)[0]]])

    # [~, ~, ~, ~, ~, ~, inside_grains, ~, ~] = grainsize_areas_planimetric(ebsd, polygon, varargin{:});

    # ASTM E-2627 takes the mean area of grains with over 100 px each, and
    # requires the average grain prior to thresholding has at least 500 px.
    # Here we allow an arbitrary number of pixel threshold.
    # px_area = polyarea(ebsd.unitCell(:,1), ebsd.unitCell(:,2));
    unitcellx = loadmat("GS_Meas\\myEBSD_high_res_1_unit_cell_x.mat")
    unitcellx = unitcellx['ebsd_unit_cell_x']
    unitcelly = loadmat("GS_Meas\\myEBSD_high_res_1_unit_cell_y.mat")
    unitcelly = unitcelly['ebsd_unit_cell_y']

    px_area = get_polygon_area_shoelace_formula(unitcellx, unitcelly)
    threshold = min_px_per_grain * px_area

    # Number of pixels per grain before threshold
    # avg_px_per_grain_before_threshold = mean(inside_grains.area / px_area);
    inside_grains_area = loadmat("GS_Meas\\myEBSD_high_res_1_inside_grains_area_before_redux.mat")
    inside_grains_area = inside_grains_area['inside_grains_area']
    avg_px_per_grain_before_threshold = np.mean(inside_grains_area / px_area)

    # Remove grains with fewer pixels than the threshold
    # excluded_grains = inside_grains(inside_grains.area >= threshold);
    # inside_grains(inside_grains.area < threshold) = [];
    # areas = inside_grains.area;
    if measurement_type == 'AsWritten':
         areas = loadmat("GS_Meas\\myEBSD_high_res_1_inside_grains_area_after_redux.mat")
         areas = areas['inside_grains_area_1']
    else:                                           # CustomMinGS
        areas = loadmat("GS_Meas\\myEBSD_high_res_1_customGS_inside_grains_area_after_redux.mat")
        areas = areas['areas']
    n = len(areas)
    Abar = np.mean(areas)

    # Calculate the number of analyzed grains per unit area
    # NOTE THAT THIS REFLECTS THE NUMBER OF GRAINS ELIMINATED BY THE THRESHOLD
    # This value is potentially useful for assessing differences between
    # E112 planimetric measurements and the E2627 standard
    # excluded_grains_area = loadmat("GS_Meas\\myEBSD_high_res_1_excluded_grains_area.mat")      # AsWritten
    # excluded_grains_area = excluded_grains_area['excluded_grains_area']                        # AsWritten
    excluded_grains_area = loadmat("GS_Meas\\myEBSD_high_res_1_customGS_excluded_grains_area.mat")      # CustomGS
    excluded_grains_area = excluded_grains_area['customGS_excluded']                                    # CustomGS
    analyzed_area = get_polygon_area_shoelace_formula(polygon[:, 0], polygon[:, 1]) - sum(excluded_grains_area)
    N_A_measured = n / analyzed_area

    G_A = G_meanbarA(Abar)

    # TODO: Plot
    # plotting subfunction
    # if ismember('PlotResults',varargin)
    #     plot(inside_grains.boundary); hold on
    #     plot(excluded_grains, 'FaceColor', 'k', 'FaceAlpha', 1.0); hold on
    #     plot(inside_grains)

    return G_A, Abar, n, N_A_measured[0], avg_px_per_grain_before_threshold, areas


def GrainSize_E112_HeynRandomLineMLI(ebsd):
    # function [G_L, lbar, n, intercept_lengths] = GrainSize_E112_HeynRandomLineMLI(ebsd, varargin)

    # [G_L, ~, ~, ~, ~, intercept_lengths, ~, ~, ~, ~, ~] = grainsize_linint_random(ebsd, 50, varargin{:});
    G_L, _, _, _, _, intercept_lengths, _, _, _, _, _ = grainsize_linint_random(ebsd)

    lbar = np.mean(intercept_lengths)
    n = len(intercept_lengths)

    return G_L, lbar, n, intercept_lengths


def GrainSize_E112_HeynRandomLinePL(ebsd):
    # [~, G_PL, ~, MIC, ~, ~, ~, ~, ~, nlines, total_line_length] = grainsize_linint_random(ebsd, 50, varargin{:});
    _, G_PL, _, MIC, _, _, _, _, _, nlines, total_line_length = grainsize_linint_random(ebsd)

    # intersection_count = total_line_length / MIC;
    intersection_count = total_line_length / MIC

    return G_PL, MIC, intersection_count, nlines, total_line_length


def GrainSize_E112_Hilliard(ebsd):
    # function [G_PL, hilliardIntCount, hilliard_lbar, circumference] = GrainSize_E112_Hilliard(ebsd, varargin)
    # %-- Hilliard Single-Circle Procedure: Count of grains intercepting a circle.
    # % Diameter should never be smaller than the largest observed grains. Test
    # % circle should be at least 3x the length of the mean lineal intercept.
    # % Recommended: test conditions that produce ~35 counts per circle.
    #
    # %-- This function generates a single circle, counts the grains intersecting the
    # % circle (E112_HilliardGrainCount.m), calculates the mean lineal intercept,
    # % and plots the data.
    # %%
    # % populate grain data
    # [grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle', 5*degree,'unitcell');
    #
    # if ismember('exclude_twins',varargin)
    #     grains = exclude_twins(grains);
    # end

    unitCell = loadmat("GS_Meas\\myEBSD_high_res_1_unit_cell.mat")
    # stepsize = 2*abs(ebsd.unitCell(1,1));
    stepsize = 2 * abs(unitCell['myUnitCell'][0][0])

    # Extract triple points
    # tP = grains.triplePoints;
    x_tP = loadmat("GS_Meas\\myEBSD_high_res_1_x_tp.mat")
    x_tP = x_tP['x_tP']
    y_tP = loadmat("GS_Meas\\myEBSD_high_res_1_y_tp.mat")
    y_tP = y_tP['y_tP']
    tpoint = np.hstack((x_tP, y_tP))
    ntpoints = np.shape(tpoint)
    ntpoints = ntpoints[0]

    # Generating the circle
    # Approximate a circle as a polygon
    # offset = 0.02; % 2pct inset from edges
    offset = 0.02  # 2 percent inset from edges
    # xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x));
    xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x))
    # ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y));
    ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y))
    # thetas = 0.0:pi/100.0:2.0*pi;
    thetas = np.arange(start=0, step=(np.pi/100), stop=((2 * np.pi) + (np.pi/100)))
    # xres = 2.0 * xcenter / length(ebsd.x);
    xres = (2 * xcenter) / len(ebsd.x)
    # yres = 2.0 * ycenter / length(ebsd.y);
    yres = (2 * ycenter) / len(ebsd.y)
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
    # hilliardPolygon = [circ_x' circ_y'];
    hilliardPolygon = np.hstack((circ_x, circ_y))

    # Generating the points where the line and the grain boundary intersect
    # hilliard_intersections = [];
    hilliard_intersections = []
    # # for n = 1:length(thetas)-1
    # for i in range(0, len(thetas)-1):
    # #     x_start = hilliardPolygon(n,1);
    #     x_start = hilliardPolygon[i][0]
    # #     x_end = hilliardPolygon(n+1,1);
    #     x_end = hilliardPolygon[i+1][0]
    # #     y_start = hilliardPolygon(n,2);
    #     y_start = hilliardPolygon[i][1]
    # #     y_end = hilliardPolygon(n+1,2);
    #     y_end = hilliardPolygon[i+1][1]
    # #     xy1 = [x_start, y_start];
    #     xy1 = [x_start, y_start]
    # #     xy2 = [x_end, y_end];
    #     xy2 = [x_end, y_end]
    #     # TODO: translate below (skipping because of MTEX use in loop)
    # #     [xi, yi] = grains.boundary.intersect(xy1, xy2);
    # #     x1 = xi(~isnan(xi));
    # #     y1 = yi(~isnan(yi));
    # #     intersect_coords = [x1',y1'];
    # #     num_int(n) = numel(x1);
    # #     hilliard_intersections = cat(1, hilliard_intersections, intersect_coords);

    num_int = loadmat("GS_Meas\\hilliard_num_int.mat")
    num_int = num_int["num_int"][0]
    hilliardIntCount = sum(num_int) - 1

    hilliard_intersections = loadmat("GS_Meas\\hilliard_intersections.mat")
    hilliard_intersections = hilliard_intersections["hilliard_intersections"]

    # Calculate the distance between intersection points and triple points
    # triplept_intersection_coordinates = [];
    triplept_intersection_coordinates = []
    # tp_thresh = 1.0; % multiples of step size
    tp_thresh = 1.0  # multiples of step size
    # for m = 1:ntpoints
    #     % distance in microns:
    #     dist = sqrt((tpoint(m,1) - hilliard_intersections(1:end,1)).^2 + ...
    #                 (tpoint(m,2) - hilliard_intersections(1:end,2)).^2) * ...
    #                  tp_thresh * stepsize;
    #
    #     %find the distance under threshold and use that as an index into xyints:
    #     coord = hilliard_intersections(dist<stepsize, :);
    #     xcoord = coord(:, 1);
    #     ycoord = coord(:, 2);
    #     triplept_intersection_coordinates = cat(1, triplept_intersection_coordinates, [xcoord, ycoord]);
    # end % triple point distance loop

    for m in range(ntpoints):
        # Distance in microns
        dist = np.sqrt(((tpoint[m, 0] - hilliard_intersections[:, 0]) ** 2) +
                       ((tpoint[m, 1] - hilliard_intersections[:, 1]) ** 2)) * tp_thresh * stepsize

        # Find the distance under threshold and use that as an index into xyints:
        current_coord = hilliard_intersections[dist < stepsize]
        if len(current_coord) != 0:
            for i in range(len(current_coord)):
                coord = current_coord[i]
                xcoord = coord[:][0]
                ycoord = coord[:][1]
                triplept_intersection_coordinates.append(np.hstack((xcoord, ycoord)))

    # Get the count of intersections through the triple points (from xcoord and ycoord)
    # xc = triplept_intersection_coordinates(:,1);
    xc = [triplept_intersection_coordinates[i][0] for i in range(len(triplept_intersection_coordinates))]
    # yc = triplept_intersection_coordinates(:,2);
    yc = [triplept_intersection_coordinates[i][1] for i in range(len(triplept_intersection_coordinates))]
    # hilliardTPcount = numel(xc)-1;
    hilliardTPcount = len(xc) - 1

    # Get the count of intersections through the triple points (from xcoord and ycoord)
    # % tpcount = numel(xcoord);

    # Add 0.5 counts for each time the line goes through a triple point
    # hilliardIntCount = hilliardIntCount + 0.5*hilliardTPcount;
    hilliardIntCount = hilliardIntCount + (0.5 * hilliardTPcount)

    # mean lineal intercept = circumference of circle / number of grains intersecting circle
    # hilliard_lbar = (2*pi*radius) / hilliardIntCount;
    hilliard_lbar = (2 * np.pi * radius) / hilliardIntCount

    # G_PL = G_meanintl(hilliardIntCount);
    G_PL = G_meanintl(hilliardIntCount)

    # circumference = 2.0 * pi * radius;
    circumference = 2.0 * np.pi * radius

    # % Plotting Subfunction
    #     if ismember('PlotResults',varargin)
    #         %--- plot the grains and grain boundaries
    #         plot(ebsd, ebsd.orientations); hold on
    #         plot(grains.boundary,'LineWidth',1); hold on
    #         %--- plot the triple points
    #         plot(tP,'color','b','linewidth',2); hold on
    #         %--- plot the circle
    #         plot(hilliardPolygon(:,1), hilliardPolygon(:,2), ...
    #             'k', 'linewidth', 3); hold on
    #         %--- plotting the intersections
    #         for n = 1:length(thetas)-1
    #             x_start = hilliardPolygon(n,1);
    #             x_end = hilliardPolygon(n+1,1);
    #             y_start = hilliardPolygon(n,2);
    #             y_end = hilliardPolygon(n+1,2);
    #             xy1 = [x_start, y_start];
    #             xy2 = [x_end, y_end];
    #             [xi, yi] = grains.boundary.intersect(xy1, xy2);
    #             x1 = xi(~isnan(xi));
    #             y1 = yi(~isnan(yi));
    #             intersect_coords = [x1',y1'];
    #             scatter(intersect_coords(:,1), ...
    #                 intersect_coords(:,2),'w','linewidth',2); hold on
    #         end % line intersection loop
    #         scatter(xc,yc,'r','linewidth',2)
    #
    #     end % varargin parse
    #
    # end % function

    return G_PL, hilliardIntCount, hilliard_lbar, circumference


def GrainSize_E112_Abrams(ebsd):
    # function [G_PL, abramsIntCount, abrams_lbar, abramsCircumference_tot] = GrainSize_E112_Abrams(ebsd)
    # Abrams Three-Circle Procedure: Three concentric and equally spaced
    # circles.
    # Placement of the three-circle test grid should yield 40-100 intersection
    # counts per slice tested.
    # If the circle intersects a triple point, the count is 2.
    # The ratio of circumference is 3:2:1

    # Populate the grains
    # TODO: Unit cell is same for this file, but might differ due to 2*degree below opposed to 5*degree (Determine diff)
    # [grains, ebsd.grainId] = calcGrains(ebsd('indexed'), 'angle', 2*degree, 'unitcell');
    # % plot(ebsd, ebsd.orientations); % plot the grains
    # % hold on
    # % plot(grains.boundary,'LineWidth',1); % plot the grain boundaries
    # % hold on

    unitCell = loadmat("GS_Meas\\myEBSD_high_res_1_unit_cell.mat")
    # stepsize = 2*abs(ebsd.unitCell(1,1));
    stepsize = 2 * abs(unitCell['myUnitCell'][0][0])  # Calculate step size for triple point calculations

    # Extract triple points
    # tP = grains.triplePoints;
    x_tP = loadmat("GS_Meas\\myEBSD_high_res_1_x_tp.mat")
    x_tP = x_tP['x_tP']
    y_tP = loadmat("GS_Meas\\myEBSD_high_res_1_y_tp.mat")
    y_tP = y_tP['y_tP']
    tpoint = np.hstack((x_tP, y_tP))
    ntpoints = np.shape(tpoint)
    ntpoints = ntpoints[0]

    # Plot triple points
    # % plot(tP,'color','b','linewidth',2); hold on

    # Plotting the largest circle:
    # Approximate a circle as a polygon
    #         offset = 0.02; % 2pct inset from edges
    offset = 0.02  # 2 percent inset from edges
    #         xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x));
    xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x))
    #         ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y));
    ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y))
    #         thetas = 0.0:pi/100.0:2.0*pi;
    thetas = np.arange(start=0, step=(np.pi / 100), stop=((2 * np.pi) + (np.pi / 100)))
    #         xres = 2.0 * xcenter / length(ebsd.x);
    xres = (2 * xcenter) / len(ebsd.x)
    #         yres = 2.0 * ycenter / length(ebsd.y);
    yres = (2 * ycenter) / len(ebsd.y)
    #         radius_lg = 0.5 * min(max(ebsd.x) - min(ebsd.x), ...
    #                        max(ebsd.y) - min(ebsd.y));
    radius_lg = 0.5 * min((max(ebsd.x) - min(ebsd.x)), (max(ebsd.y) - min(ebsd.y)))
    #         inset = max(numel(ebsd.x) * offset * xres, numel(ebsd.y) * offset * yres);
    inset = max((len(ebsd.x) * offset * xres), (len(ebsd.y) * offset * yres))
    #         radius_lg = radius_lg - inset; % inset from the edges of the scan
    radius_lg = radius_lg - inset  # Inset from the edges of the scan
    #         circumference_lg = 2 * pi * radius_lg;
    circumference_lg = 2 * np.pi * radius_lg
    #         circ_x_lg = radius_lg * cos(thetas) + xcenter;
    circ_x_lg = radius_lg * np.cos(thetas) + xcenter
    #         circ_y_lg = radius_lg * sin(thetas) + ycenter;
    circ_y_lg = radius_lg * np.sin(thetas) + ycenter
    #         polygon_lg = [circ_x_lg' circ_y_lg']; % x/y coords of each line segment
    polygon_lg = np.hstack((circ_x_lg, circ_y_lg))  # x/y coords of each line segment

    #         % plot the largest circle
    # %         plot(circ_x_lg,circ_y_lg, 'k', 'linewidth', 3)
    # %         hold on

    # Extract the grain boundary/circle intersection data
    #         g = size(polygon_lg);
    g = np.shape(polygon_lg)
    #         gg = g(1);
    gg = g[0]

    #         abrams_intersections_lg = [];
    abrams_intersections_lg = []
    # #         for n = 1:gg-1
    # for i in range(gg - 1):
    # #             x_start = polygon_lg(n,1);
    #     x_start = polygon_lg[i][0]
    # #             x_end = polygon_lg(n+1,1);
    #     x_end = polygon_lg[i + 1][0]
    # #             y_start = polygon_lg(n,2);
    #     y_start = polygon_lg[i][1]
    # #             y_end = polygon_lg(n+1,2);
    #     y_end = polygon_lg[i + 1][1]
    # #             xy1 = [x_start, y_start];
    #     xy1 = [x_start, y_start]
    # #             xy2 = [x_end, y_end];
    #     xy2 = [x_end, y_end]
    #               # TODO: translate below (skipping because of MTEX use in loop)
    # #             [xi, yi] = grains.boundary.intersect(xy1, xy2);
    # #             x1 = xi(~isnan(xi));
    # #             y1 = yi(~isnan(yi));
    # #             intersect_coords_lg = [x1',y1'];
    # #             %scatter(intersect_coords_lg(:,1),intersect_coords_lg(:,2),'w','linewidth',2)
    # #             num_int_lg(n) = numel(x1);
    # #             abrams_intersections_lg = cat(1, abrams_intersections_lg, intersect_coords_lg);
    # #         end % line intersection loop
    # #         abramsIntCount_lg = sum(num_int_lg)-1;

    abrams_intersections_lg = loadmat("GS_Meas\\abrams_intersections_lg.mat")
    abrams_intersections_lg = abrams_intersections_lg["abrams_intersections_lg"]
    num_int_lg = loadmat("GS_Meas\\abrams_num_int_lg.mat")
    num_int_lg = num_int_lg["num_int_lg"][0]
    abramsIntCount_lg = sum(num_int_lg) - 1

    # Plotting the medium circle:
    # Approximate a circle as a polygon
    #         % offset = 0.02; % 2pct inset from edges
    #         xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x));
    #         ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y));
    #         thetas = 0.0:pi/100.0:2.0*pi;
    #         % xres = 2.0 * xcenter / length(ebsd.x);
    #         % yres = 2.0 * ycenter / length(ebsd.y);
    #         circumference_med = circumference_lg / 1.5;
    circumference_med = circumference_lg / 1.5
    #         radius_med = circumference_med / (2 * pi);
    radius_med = circumference_med / (2 * np.pi)
    #         % inset = max(numel(ebsd.x) * offset * xres, numel(ebsd.y) * offset * yres);
    #         radius_med = radius_med - inset; % inset from the edges of the scan
    radius_med = radius_med - inset  # Inset from the edges of the scan
    #         circ_x_med = radius_med * cos(thetas) + xcenter;
    circ_x_med = radius_med * np.cos(thetas) + xcenter
    #         circ_y_med = radius_med * sin(thetas) + ycenter;
    circ_y_med = radius_med * np.sin(thetas) + ycenter
    #         polygon_med = [circ_x_med' circ_y_med'];
    polygon_med = np.hstack((circ_x_med, circ_y_med))

    # %         plot(circ_x_med,circ_y_med, 'k', 'linewidth', 3)
    # %         hold on

    #         g = size(polygon_med);
    g = np.shape(polygon_med)
    #         gg = g(1);
    gg = g[0]

    #         abrams_intersections_med = [];
    abrams_intersections_med = []
    # #         for n = 1:gg-1
    # for i in range(gg - 1):
    # #             x_start = polygon_med(n,1);
    #     x_start = polygon_med[i][0]
    # #             x_end = polygon_med(n+1,1);
    #     x_end = polygon_med[i + 1][0]
    # #             y_start = polygon_med(n,2);
    #     y_start = polygon_med[i][1]
    # #             y_end = polygon_med(n+1,2);
    #     y_end = polygon_med[i + 1][1]
    # #             xy1 = [x_start, y_start];
    #     xy1 = [x_start, y_start]
    # #             xy2 = [x_end, y_end];
    #     xy2 = [x_end, y_end]
    #               # TODO: translate below (skipping because of MTEX use in loop)
    # #             [xi, yi] = grains.boundary.intersect(xy1, xy2);
    # #             x1 = xi(~isnan(xi));
    # #             y1 = yi(~isnan(yi));
    # #             intersect_coords_med = [x1',y1'];
    # #             %scatter(intersect_coords_med(:,1),intersect_coords_med(:,2),'w','linewidth',2)
    # #             num_int_med(n) = numel(x1);
    # #             abrams_intersections_med = cat(1, abrams_intersections_med, intersect_coords_med);
    # #         end % line intersection loop
    # #         abramsIntCount_med = sum(num_int_med)-1;

    abrams_intersections_med = loadmat("GS_Meas\\abrams_intersections_med.mat")
    abrams_intersections_med = abrams_intersections_med["abrams_intersections_med"]
    num_int_med = loadmat("GS_Meas\\abrams_num_int_med.mat")
    num_int_med = num_int_med["num_int_med"][0]
    abramsIntCount_med = sum(num_int_med) - 1

    # Plotting the smallest circle
    # Approximate a circle as a polygon
    #         % offset = 0.02; % 2pct inset from edges
    #         xcenter = 0.5 * (max(ebsd.x) - min(ebsd.x));
    #         ycenter = 0.5 * (max(ebsd.y) - min(ebsd.y));
    #         thetas = 0.0:pi/100.0:2.0*pi;
    #         % xres = 2.0 * xcenter / length(ebsd.x);
    #         % yres = 2.0 * ycenter / length(ebsd.y);
    #         circumference_sm = circumference_lg / 3;
    circumference_sm = circumference_lg / 3
    #         radius_sm = circumference_sm / (2 * pi);
    radius_sm = circumference_sm / (2 * np.pi)
    #         % inset = max(numel(ebsd.x) * offset * xres, numel(ebsd.y) * offset * yres);
    #         radius_sm = radius_sm - inset; % inset from the edges of the scan
    radius_sm = radius_sm - inset  # Inset from the edges of the scan
    #         circ_x_sm = radius_sm * cos(thetas) + xcenter;
    circ_x_sm = radius_sm * np.cos(thetas) + xcenter
    #         circ_y_sm = radius_sm * sin(thetas) + ycenter;
    circ_y_sm = radius_sm * np.sin(thetas) + ycenter
    #         polygon_sm = [circ_x_sm' circ_y_sm'];
    polygon_sm = np.hstack((circ_x_sm, circ_y_sm))

    # %         plot(circ_x_sm,circ_y_sm, 'k', 'linewidth', 3)
    # %         hold on
    #
    #         g = size(polygon_sm);
    g = np.shape(polygon_sm)
    #         gg = g(1);
    gg = g[0]

    #         abrams_intersections_sm = [];
    abrams_intersections_sm = []
    # #         for n = 1:gg-1
    # for i in range(gg - 1):
    # #             x_start = polygon_sm(n,1);
    # x_start = polygon_sm[i][0]
    # #             x_end = polygon_sm(n+1,1);
    # x_end = polygon_sm[i + 1][0]
    # #             y_start = polygon_sm(n,2);
    # y_start = polygon_sm[i][1]
    # #             y_end = polygon_sm(n+1,2);
    # y_end = polygon_sm[i + 1][1]
    # #             xy1 = [x_start, y_start];
    # xy1 = [x_start, y_start]
    # #             xy2 = [x_end, y_end];
    # xy2 = [x_end, y_end]
    #                # TODO: translate below (skipping because of MTEX use in loop)
    # #             [xi, yi] = grains.boundary.intersect(xy1, xy2);
    # #             x1 = xi(~isnan(xi));
    # #             y1 = yi(~isnan(yi));
    # #             intersect_coords_sm = [x1',y1'];
    # # %             scatter(intersect_coords_sm(:,1),intersect_coords_sm(:,2),'w','linewidth',2)
    # #             num_int_sm(n) = numel(x1);
    # #             abrams_intersections_sm = cat(1, abrams_intersections_sm, intersect_coords_sm);
    # #         end % line intersection loop
    # #         abramsIntCount_sm = sum(num_int_sm)-1;

    abrams_intersections_sm = loadmat("GS_Meas\\abrams_intersections_sm.mat")
    abrams_intersections_sm = abrams_intersections_sm["abrams_intersections_sm"]
    num_int_sm = loadmat("GS_Meas\\abrams_num_int_sm.mat")
    num_int_sm = num_int_sm["num_int_sm"][0]
    abramsIntCount_sm = sum(num_int_sm) - 1

    # Concatenating polygon data
    # % abramsPolygon = cat(1, abramsPolygon, polygon_sm, polygon_med, polygon_lg);
    # % g = size(abramsPolygon);
    # % gg = g(1);

    # abrams_intersections = cat(1, abrams_intersections_sm, abrams_intersections_med, ...
    #     abrams_intersections_lg);

    abrams_intersections = np.vstack([abrams_intersections_sm, abrams_intersections_med, abrams_intersections_lg])

    # Calculate the distance between intersection points and triple points
    #     triplept_intersection_coordinates = [];
    triplept_intersection_coordinates = []
    #     tp_thresh = 1.0; % multiples of step size
    tp_thresh = 1.0  # Multiples of step size
    #     for m = 1:ntpoints
    #         % distance in microns:
    #         dist = sqrt((tpoint(m,1) - abrams_intersections(1:end,1)).^2 + ...
    #                     (tpoint(m,2) - abrams_intersections(1:end,2)).^2) * ...
    #                      tp_thresh * stepsize;
    #
    #         %find the distance under threshold and use that as an index into xyints:
    #         coord = abrams_intersections(dist<stepsize, :);
    #         xcoord = coord(:, 1);
    #         ycoord = coord(:, 2);
    #         triplept_intersection_coordinates = cat(1, triplept_intersection_coordinates, [xcoord, ycoord]);
    #     end % triple point distance loop

    for m in range(ntpoints):
        # Distance in microns
        dist = np.sqrt(((tpoint[m, 0] - abrams_intersections[:, 0]) ** 2) +
                       ((tpoint[m, 1] - abrams_intersections[:, 1]) ** 2)) * tp_thresh * stepsize

        # Find the distance under threshold and use that as an index into xyints:
        current_coord = abrams_intersections[dist < stepsize]
        if len(current_coord) != 0:
            for i in range(len(current_coord)):
                coord = current_coord[i]
                xcoord = coord[:][0]
                ycoord = coord[:][1]
                triplept_intersection_coordinates.append(np.hstack((xcoord, ycoord)))

    # Get the count of intersections through the triple points (from xcoord and ycoord)
    # xc = triplept_intersection_coordinates(:,1);
    xc = [triplept_intersection_coordinates[i][0] for i in range(len(triplept_intersection_coordinates))]
    # yc = triplept_intersection_coordinates(:,2);
    yc = [triplept_intersection_coordinates[i][1] for i in range(len(triplept_intersection_coordinates))]
    # % hold on
    # % scatter(xc,yc,'r','linewidth',2)

    # abramsTPcount = numel(xc)-1;
    abramsTPcount = len(xc) - 1

    # abramsIntCount = abramsIntCount_sm + abramsIntCount_med + abramsIntCount_lg;
    abramsIntCount = abramsIntCount_sm + abramsIntCount_med + abramsIntCount_lg
    # abramsIntCount = abramsIntCount + abramsTPcount;
    abramsIntCount = abramsIntCount + abramsTPcount

    # Total line length = total circumference of circles
    # abramsCircumference_tot = circumference_lg + circumference_med + circumference_sm;
    abramsCircumference_tot = circumference_lg + circumference_med + circumference_sm

    # N_L = abramsIntCount / abramsCircumference_tot;
    N_L = abramsIntCount / abramsCircumference_tot
    # abrams_lbar = 1/N_L;
    abrams_lbar = 1/N_L

    # G_PL = G_meanintl(N_L);
    G_PL = G_meanintl(N_L)

    return G_PL, abramsIntCount, abrams_lbar, abramsCircumference_tot


def GrainSize_E930_ALA(ebsd, G2):
    # function [G_largestGrain, volFraction] = GrainSize_E930_ALA(ebsd, G2, varargin)
    # ASTM Standard E930 Procedure 6.4 (Referee Procedure for Image Analysis)

    # G2 ---G number calculated from G_meanbarA function

    # [grains, ebsd.grainId] = calcGrains(ebsd('indexed'),'angle', 5*degree,'unitcell');
    unitCell = loadmat("GS_Meas\\myEBSD_high_res_1_unit_cell.mat")

    # if ismember('exclude_twins',varargin)
    #     grains = exclude_twins(grains);
    # end

    # areas = grains.area;
    areas = loadmat("GS_Meas\\myEBSD_high_res_1_area.mat")
    areas = areas['areas']

    # Calculating the ASTM G number for the largest grain detected
    # extract the largest area detected
    maxGrainArea = max(areas)

    # calculate the ASTM G number for the area in square microns
    A = 2.0 * np.log2(254.0) + 1.0
    B = 1.0 / (np.log10(2.0))
    G_largestGrain = A - B * np.log10(maxGrainArea)

    # G number that represents the average grain size from average area
    G_scanAvg = G2
    G_scanAvg = G_scanAvg - 3.0  # Subtract 3 ASTM numbers

    # convert G_scanAvg back to area
    calcArea = 10.0 ** ((G_scanAvg - A) / (-1 * B))

    # Sum of the areas larger than 'calcArea'
    # largestGrains = sum(areas(areas>calcArea,:));
    large_areas = [area for area in areas if area > calcArea]
    largestGrains = sum(large_areas)

    # Calculating the volume fraction
    volFraction = largestGrains / sum(areas)
    # If largest grains make up more than 5% of the area, output a warning
    if volFraction > 0.05:
        warnings.warn(f"Volume Fraction greater than 0.05")

    return G_largestGrain, volFraction


def main():
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

    slice_index = 2
    plane_normal = 'z'
    load_d3d_slice(user_dstruct, slice_index, plane_normal)

    # Make all slices have exactly the same dimensions for comparable grain sizes
    res_adjust = 200.0 / (max(myEBSD.x) - min(myEBSD.x))
    myEBSD.x = res_adjust * myEBSD.x
    myEBSD.y = res_adjust * myEBSD.y

    # pickle_out = open(f"myHighResEBSD.pickle", "wb")
    # pickle.dump(myEBSD, pickle_out)
    # pickle_out.close()

    # Do some grain size measurements!
    # G_S, N_A_S, n_S = GrainSize_E112_SaltikovPlanimetric(myEBSD)
    # print(G_S, N_A_S, n_S)
    # G_J, N_A_J, n_J = GrainSize_E112_JeffriesPlanimetric(myEBSD)
    # print(G_J, N_A_J, n_J)
    # G_A1, Abar_A1, n_A1, N_A_measured_A1, avg_px_per_grain_after_threshold, areas_A1 = \
    #     GrainSize_E2627_AsWritten(myEBSD)
    # print(G_A1, Abar_A1, n_A1, N_A_measured_A1, avg_px_per_grain_after_threshold, areas_A1)
    # G_A2, Abar_A2, n_A2, N_A_measured_A2, avg_px_per_grain_before_threshold, areas_A2 = \
    #     GrainSize_E2627_CustomMinGS(myEBSD)
    # print(G_A2, Abar_A2, n_A2, N_A_measured_A2, avg_px_per_grain_before_threshold, areas_A2)
    # TODO: Incomplete translation due to randlin function
    G_L, lbar, n_L_intercepts, intercept_lengths_L = GrainSize_E112_HeynRandomLineMLI(myEBSD)
    print(G_L, lbar, n_L_intercepts, intercept_lengths_L)
    # TODO: Incomplete translation due to randlin function
    # G_PL, P_L, PL_intersection_count, nlines, Heyn_total_line_length = GrainSize_E112_HeynRandomLinePL(myEBSD)
    # print(G_PL, P_L, PL_intersection_count, nlines, Heyn_total_line_length)
    # TODO: Incomplete translation due to MTEX interaction in for loop
    # G_Hilliard, hilliardIntCount, hilliard_lbar, hilliardCircumference = GrainSize_E112_Hilliard(myEBSD)
    # print(f"G_Hilliard = {G_Hilliard}, hilliardIntCount = {hilliardIntCount}, hilliard_lbar = {hilliard_lbar}, "
    #       f"hilliardCircumference = {hilliardCircumference}")
    # TODO: Incomplete translation due to MTEX interaction in for loop
    # G_Abrams, abramsIntCount, abrams_lbar, abramsCircumference = GrainSize_E112_Abrams(myEBSD)
    # print(f"G_Abrams = {G_Abrams}, abramsIntCount = {abramsIntCount}, abrams_lbar = {abrams_lbar}, "
    #       f"abramsCircumference = {abramsCircumference}")
    # G_largestGrain, volFraction = GrainSize_E930_ALA(myEBSD, G_S)
    # print(G_largestGrain, volFraction)


if __name__ == '__main__':
    main()
