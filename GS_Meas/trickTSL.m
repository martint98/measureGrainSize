function trickTSL(ebsd, name)
ttt = tic;
% this script takes in the orientations a list of row vectors of
% bunge angles and then generates a fake
% ANG file to TSL OIM v5.3 format specs.
% run this script in a directory containing a list of Bunge angles in a directory
% and it will save a *.ang file with HCP symmetry that 
%
% A.L. Pilchak, D. Banerjee, J.C. Williams, Observations of grain boundary alpha and alpha sideplate formation in alpha + beta titanium
% pilchak@matsceng.ohio-state.edu
%  
%
%% find the file names
%
% start a clock
%tic
symm = 'HCP'; % a string, options: 'FCC', 'BCC', 'HCP'
%fnames = dir('*.TEX2');
%
%[numInFiles dummy1] = size(fnames);
%phase = 0; % 
%
%
%% import the data, one file at a time
%
%for jj = 1:numInFiles
%
%
%    [path,name,ext] = fileparts(fnames(jj).name);
%    
%    fname = sprintf('%s',name,ext);
%
    ext2 = '.ang';
    WriteFname = sprintf('%s',name, ext2);
%
%
%  update the user
%    fprintf('\r\rWorking on file number (%s) of (%s), %s \r',num2str(jj),num2str(numInFiles),fname)
%    
%
% You could also input the files by modifying the code and having the user enter a list of ROW vector Euler angles [phi1 PHI phi2;];
%    BungeAngles = load(fname);
%    

% Preallocate variables that need to be obtained from ebsd object
N = size(ebsd, 1);
BungeAngles = zeros(N, 3);
Xpositions = zeros(N, 1);
Ypositions = zeros(N, 1);
phase = zeros(N, 1);

phase(1:N, 1) = ebsd(1:N).phase + 1;
for i=1:N % This loop is VERY inefficient but necessary due to the way mtex stores data in the ebsd object
    if phase(i) > 0
        BungeAngles(i, 1) = ebsd(i).orientations.phi1;
        BungeAngles(i, 2) = ebsd(i).orientations.Phi;
        BungeAngles(i, 3) = ebsd(i).orientations.phi2;
    else
        BungeAngles(i, 1) = 0.0;
        BungeAngles(i, 2) = 0.0;
        BungeAngles(i, 3) = 0.0;        
    end
end
Xpositions(1:N, 1) = ebsd(1:N).prop.x;
Ypositions(1:N, 1) = ebsd(1:N).prop.y;

    
%
% convert to radians
bunge_radians = pi/180 .* BungeAngles(:,1:3);
% bunge_radians =  BungeAngles(:,1:3);
%
% make some fake x,y coords
%Xpositions = [1:1:N]';
%Ypositions = [1:1:N]';
%
% generate the matrix of space fillers
filler = ones(N,1);
%
% create the matrix containing all elements required for TSL file
% the orientations are real and the phase is real (1) for Ti-Alpha
%
TSLfile = [bunge_radians Xpositions Ypositions filler(:,1) filler(:,1) phase.*filler(:,1) filler(:,1) filler(:,1)];

fileWriteID = fopen(WriteFname,'w+');

if strcmp(symm,'FCC') == 1
            %
            fprintf(fileWriteID,'# TEM_PIXperUM          1.000000\r\n');
            fprintf(fileWriteID,'# x-star                0.571337\r\n');
            fprintf(fileWriteID,'# y-star                0.882437\r\n');
            fprintf(fileWriteID,'# z-star                0.657823\r\n');
            fprintf(fileWriteID,'# WorkingDistance       19.000000\r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# Phase 1\r\n');
            fprintf(fileWriteID,'# MaterialName  	Nickel\r\n');
            fprintf(fileWriteID,'# Formula     	Ni\r\n');
            fprintf(fileWriteID,'# Info 		\r\n');
            fprintf(fileWriteID,'# Symmetry              43\r\n');
            fprintf(fileWriteID,'# LatticeConstants      3.560 3.560 3.560  90.000  90.000  90.000\r\n');
            fprintf(fileWriteID,'# NumberFamilies        69\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -1 -1 1 11.936152 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -2  0 1 10.514876 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -2  2 1 7.489179 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -3 -1 1 6.273238 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -2 -2 0 5.956645 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -4  0 0 4.956584 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -3 -3 0 4.425857 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -4  2 0 4.271905 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -4 -2 0 3.737260 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -5 -1 0 3.423848 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -3 -3 0 3.423848 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -4  4 0 2.963146 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -5 -3 0 2.768213 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -4 -4 0 2.708809 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -6  0 0 2.708809 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -6  2 0 2.479068 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -5 -3 0 2.335206 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -6 -2 0 2.292190 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -4 -4 0 2.124829 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -5 -5 0 2.006490 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -7 -1 0 2.006490 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -6  4 0 1.975567 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -6 -4 0 1.854753 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -7 -3 0 1.766952 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -5 -5 0 1.766952 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -8  0 0 1.642427 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -7 -3 0 1.578907 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -8  2 0 1.558051 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -6 -4 0 1.558051 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -6  6 0 1.476123 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -8 -2 0 1.476123 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -5 -5 0 1.423998 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -7 -5 0 1.423998 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -6 -6 0 1.408157 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -8  4 0 1.345813 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -9 -1 0 1.300071 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -7 -5 0 1.300071 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -8 -4 0 1.285009 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -6 -6 0 1.231948 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -9 -3 0 1.196316 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -8 -4 0 1.138211 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -9 -3 0 1.104073 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -7 -7 0 1.104073 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -7 -5 0 1.104073 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -8  6 0 1.094247 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -8 -6 0 1.058645 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -7 -7 0 1.032390 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -9 -5 0 1.032390 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 6 -6 -6 0 1.023721 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -9 -5 0 0.965601 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -8 -6 0 0.958806 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -7 -7 0 0.912034 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -8  8 0 0.879436 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -9 -5 0 0.861084 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -9 -7 0 0.861084 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -8 -8 0 0.855450 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 6 -8 -6 0 0.833125 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -9 -7 0 0.816595 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -8 -8 0 0.789438 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 7 -7 -7 0 0.773942 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -9 -7 0 0.741021 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -9 -9 0 0.708940 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 6 -8 -8 0 0.704986 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -9 -9 0 0.680567 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 7 -9 -7 0 0.653385 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -9 -9 0 0.628454 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 8 -8 -8 0 0.614170 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 7 -9 -9 0 0.561518 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 9 -9 -9 0 0.477942 0\r\n');
            fprintf(fileWriteID,'# Categories0 0 0 0 0 \r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# GRID: SqrGrid\r\n');
            fprintf(fileWriteID,'# XSTEP: 0.650000\r\n');
            fprintf(fileWriteID,'# YSTEP: 0.650000\r\n');
            fprintf(fileWriteID,'# NCOLS_ODD: 55\r\n');
            fprintf(fileWriteID,'# NCOLS_EVEN: 55\r\n');
            fprintf(fileWriteID,'# NROWS: 48\r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# OPERATOR: 	AnyStitch\r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# SAMPLEID: \r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# SCANID: 	\r\n');
            fprintf(fileWriteID,'#\r\n');
            %    
       elseif strcmp(symm,'HCP') == 1
            %
            fprintf(fileWriteID,'# TEM_PIXperUM          1.000000\r\n');
            fprintf(fileWriteID,'# x-star                0.568837\r\n');
            fprintf(fileWriteID,'# y-star                0.873537\r\n');
            fprintf(fileWriteID,'# z-star                0.658923\r\n');
            fprintf(fileWriteID,'# WorkingDistance       25.000000\r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# Phase 1\r\n');
            fprintf(fileWriteID,'# MaterialName  	Titanium (Alpha)\r\n');
            fprintf(fileWriteID,'# Formula     	Ti\r\n');
            fprintf(fileWriteID,'# Info 		\r\n');
            fprintf(fileWriteID,'# Symmetry              62\r\n');
            fprintf(fileWriteID,'# LatticeConstants      2.950 2.950 4.680  90.000  90.000 120.000\r\n');
            fprintf(fileWriteID,'# NumberFamilies        8\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1  0  0 1 0.000000 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0  0  2 1 0.000000 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1  0  1 1 0.000000 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1  0  2 1 0.000000 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1  1  0 1 0.000000 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1  0  3 1 0.000000 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1  1  2 1 0.000000 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2  0  1 1 0.000000 1\r\n');
            fprintf(fileWriteID,'# Categories0 0 0 0 0 \r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# GRID: SqrGrid\r\n');
            fprintf(fileWriteID,'# XSTEP: 0.400000\r\n'); % these don't really matter as they aren't used by TSL
            fprintf(fileWriteID,'# YSTEP: 0.400000\r\n'); % although new things like "crystal overlay size" haven't been tested.
            fprintf(fileWriteID,'# NCOLS_ODD: 188\r\n');
            fprintf(fileWriteID,'# NCOLS_EVEN: 187\r\n');
            fprintf(fileWriteID,'# NROWS: 217\r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# OPERATOR: AnyStitch	\r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# SAMPLEID: \r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# SCANID: 	\r\n');
            fprintf(fileWriteID,'#\r\n');    
           %
          elseif strcmp(symm,'BCC') == 1
            fprintf(fileWriteID,'# TEM_PIXperUM          1.000000 \r\n');
            fprintf(fileWriteID,'# x-star                0.516878\r\n');
            fprintf(fileWriteID,'# y-star                0.654829\r\n');
            fprintf(fileWriteID,'# z-star                0.577533\r\n');
            fprintf(fileWriteID,'# WorkingDistance       25.000000\r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# Phase 1\r\n');
            fprintf(fileWriteID,'# MaterialName  	Titanium - Beta\r\n');
            fprintf(fileWriteID,'# Formula     	Ti\r\n');
            fprintf(fileWriteID,'# Info 		\r\n');
            fprintf(fileWriteID,'# Symmetry              43\r\n');
            fprintf(fileWriteID,'# LatticeConstants      3.240 3.240 3.240  90.000  90.000  90.000\r\n');
            fprintf(fileWriteID,'# NumberFamilies        100\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -1  1 1 6.689197 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -2  0 1 4.647529 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -2 -1 1 3.676836 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -2  2 0 3.063804 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -3  1 1 2.628516 1\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -2 -2 0 2.308226 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -3 -2 0 2.050676 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -4  0 0 1.840713 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -3  3 0 1.673577 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -4 -1 0 1.673577 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -4  2 0 1.534913 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -3 -3 0 1.416775 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -4 -2 0 1.323520 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -5  1 0 1.234076 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -4 -3 0 1.234076 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -5 -2 0 1.087964 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -4  4 0 1.030439 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -5  3 0 0.974686 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -4 -3 0 0.974686 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -6  0 0 0.924103 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -4 -4 0 0.924103 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -5 -3 0 0.885359 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -6 -1 0 0.885359 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -6  2 0 0.847621 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -5 -4 0 0.810816 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -6 -2 0 0.775638 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -6 -3 0 0.749065 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -4 -4 0 0.723063 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -7  1 0 0.697598 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -5 -4 0 0.697598 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -5  5 0 0.697598 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -6  4 0 0.672637 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -6 -3 0 0.650963 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -5 -5 0 0.650963 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -7 -2 0 0.650963 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -6 -4 0 0.632631 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -7  3 0 0.614623 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -6 -5 0 0.579516 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -7 -3 0 0.579516 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -8  0 0 0.563909 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -8 -1 0 0.550419 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -7 -4 0 0.550419 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -5 -5 0 0.550419 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -6 -4 0 0.537131 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -8  2 0 0.537131 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -6 -5 0 0.524037 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -8 -2 0 0.511129 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -6  6 0 0.511129 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -7  5 0 0.498400 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -8 -3 0 0.498400 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -7 -4 0 0.498400 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -6 -6 0 0.487867 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -7 -5 0 0.477538 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -8  4 0 0.467341 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -8 -3 0 0.457270 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -9  1 0 0.457270 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -8 -4 0 0.447322 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -7 -6 0 0.437593 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -9 -2 0 0.437593 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -6 -5 0 0.437593 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -6 -6 0 0.429819 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -8 -5 0 0.422133 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -9  3 0 0.422133 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -7 -5 0 0.422133 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -9 -3 0 0.407013 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -7 -6 0 0.407013 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -8 -4 0 0.399574 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -9 -4 0 0.392212 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -7  7 0 0.392212 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -8 -5 0 0.392212 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -8  6 0 0.385770 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -7 -7 0 0.379758 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -8 -6 0 0.373804 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -9  5 0 0.367908 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -9 -4 0 0.367908 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 6 -6 -6 0 0.362066 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -7 -6 0 0.356279 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -9 -5 0 0.356279 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -8 -7 0 0.345373 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -8 -5 0 0.345373 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -7 -7 0 0.345373 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -8 -6 0 0.340302 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -9 -6 0 0.335275 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -8 -7 0 0.325346 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -9 -5 0 0.325346 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -9 -6 0 0.315579 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -8  8 0 0.311475 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -9  7 0 0.307751 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -8 -8 0 0.304055 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 2 -9 -7 0 0.300388 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 6 -7 -7 0 0.300388 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 6 -8 -6 0 0.296747 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -8 -7 0 0.293134 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 5 -9 -6 0 0.285984 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -8 -8 0 0.282447 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 1 -9 -8 0 0.278935 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 4 -9 -7 0 0.278935 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 3 -9 -8 0 0.265119 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 0 -9  9 0 0.252565 0\r\n');
            fprintf(fileWriteID,'# hklFamilies   	 7 -8 -7 0 0.252565 0\r\n');
            fprintf(fileWriteID,'# Categories0 0 0 0 0 \r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# GRID: SqrGrid \r\n');
            fprintf(fileWriteID,'# XSTEP: 1.000000\r\n');
            fprintf(fileWriteID,'# YSTEP: 1.000000\r\n');
            fprintf(fileWriteID,'# NCOLS_ODD: 181\r\n');
            fprintf(fileWriteID,'# NCOLS_EVEN: 181\r\n');
            fprintf(fileWriteID,'# NROWS: 181\r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# OPERATOR: 	AnyStitch\r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# SAMPLEID: \r\n');
            fprintf(fileWriteID,'#\r\n');
            fprintf(fileWriteID,'# SCANID: 	\r\n');
            fprintf(fileWriteID,'#\r\n');

end 
    
fclose all;

% write the file
dlmwrite(WriteFname,TSLfile,'delimiter',' ','-append','precision','%.4f');
clear BungeAngles TSLfile 

%end
toc(ttt)