function [data, dstruct] = load_d3d_file(fpth)

%%--- Import h5 data from structure file ---%
% Read the dataset
data = h5read(fpth, '/DataContainers/SyntheticVolumeDataContainer/CellData/EulerAngles');

% Extract synthetic structure generation stats
dstruct.mu = h5read(fpth, '/DataContainers/StatsGeneratorDataContainer/CellEnsembleData/Statistics/1/FeatureSize Distribution/Average');
dstruct.sigma = h5read(fpth, '/DataContainers/StatsGeneratorDataContainer/CellEnsembleData/Statistics/1/FeatureSize Distribution/Standard Deviation');
tmp = h5read(fpth, '/DataContainers/StatsGeneratorDataContainer/CellEnsembleData/Statistics/1/BinNumber');
dstruct.numbins = length(tmp);
dstruct.minsize = min(tmp);
dstruct.maxsize = max(tmp);
dstruct.meansize = exp(dstruct.mu + 0.5*dstruct.sigma^2);
dstruct.binstepsize = mean(tmp(2:end)-tmp(1:end-1));
dstruct.res = h5read(fpth, '/DataContainers/SyntheticVolumeDataContainer/_SIMPL_GEOMETRY/SPACING');
dstruct.dims = double(h5read(fpth, '/DataContainers/SyntheticVolumeDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));
dstruct.mincutoff = (log(dstruct.minsize) - dstruct.mu) / dstruct.sigma;
dstruct.maxcutoff = (log(dstruct.maxsize) - dstruct.mu) / dstruct.sigma;
dstruct.ngrains = double(h5readatt(fpth, '/DataContainers/SyntheticVolumeDataContainer/Grain Data', 'TupleDimensions'));

end