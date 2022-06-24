function ebsd = load_d3d_slice(data, sliceindex, plane_normal)

shape = size(data);
xdim = shape(2);
ydim = shape(3);
zdim = shape(4); % UNUSED NOW, NEEDED LATER. SEE TO-DO BELOW.

% Determine which slice to be analyzed from a certain direction
if plane_normal == 'z'
    slice = squeeze(data(:, :, :, sliceindex));
elseif plane_normal == 'y'
    slice = squeeze(data(:, :, sliceindex, :));
elseif plane_normal == 'x'
    slice = squeeze(data(:, sliceindex, :, :));
end

% TO-DO: The reshaping needs to be generalized for non-cube datasets!!!

% Reshape the slice into an n x 3 array of presumed Euler angles
Eul = reshape(slice, 3, [])';
% Generate an nx1 array of x-coordinates
x = repmat(1:xdim, [1, ydim])';
% Generate an nx1 array of y-coordinates
y = reshape(repmat(1:ydim, [xdim, 1]), 1, [])';
% Generate an nx1 array of phase indices
phase = ones(xdim * ydim, 1);

%--- Import the slice into MTEX ---%
% crystal symmetry
CS = {...
    'notIndexed',...
    crystalSymmetry('m-3m', [3 3 3], 'mineral', 'genericCubic', 'color', [0 0 1]),...
    'notIndexed'};

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

q = rotation.byEuler(Eul(:,1), Eul(:,2), Eul(:,3));
prop.x = x;
prop.y = y;
ebsd = EBSD(q, phase, CS, prop);

end