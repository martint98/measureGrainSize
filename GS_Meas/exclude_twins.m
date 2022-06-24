function grains_excluding_twins = exclude_twins(grains_including_twins)
% TODO: Extend to deformation twins in BCC and HCP

% Find & merge FCC annealing twins
CS = grains_including_twins.CS;
twinning = orientation.map(Miller(1,1,-1,CS),Miller(-1,1,-1,CS),Miller(-1,1,1,CS));
gB = grains_including_twins.boundary('indexed','indexed');
isTwinning = angle(gB.misorientation,twinning) < 5*degree; %TODO! Does 5 deg make sense? Make variable?
twinBoundary = gB(isTwinning);
[grains_excluding_twins,~] = merge(grains_including_twins,twinBoundary);