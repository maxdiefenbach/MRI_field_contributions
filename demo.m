clear all; close all; clc;
addpath('utils/')
addpath(genpath('fwtoolbox_v1_code/'))  % ISMRM water-fat toolbox


%% load data

% load mat-file with imaging data after reconstruction
load('20151101_151725_0302_ImDataParams.mat')
signal0 = ImDataParams.signal;
signal = signal0;
TE_s = ImDataParams.TE_s;
centerFreq_Hz = ImDataParams.centerFreq_Hz;
B0dir = ImDataParams.B0dir;
voxelSize_mm = ImDataParams.voxelSize_mm;
isl = 37;
p0 = angle(signal(:, :, isl, end));

% load mat-file with coefficients of the spherical harmonic expantion of
% the magnet inhomogeneities and the shimfield 
load('20151101_151725_0302_B0params.mat')
B0params


%% estimate and demodulate field contributions

% magnet inhomogeneities
sz = size(signal);
matrixSize = sz(1:3);
transform = B0params.affMat_ijk2xyMz;
magnetInhomogeneities_Hz = get_magnetInhomogeneities_Hz(B0params, matrixSize, transform);
signal = demodulate_field_Hz(signal, TE_s, magnetInhomogeneities_Hz);
p1 = angle(signal(:, :, isl, end));

% shim field
shimField_Hz = get_shimField_Hz(B0params.shimValues, matrixSize, transform);
signal = demodulate_field_Hz(signal, TE_s, shimField_Hz);
p2 = angle(signal(:, :, isl, end));

% object-based fieldmap estimate
objectBasedFieldmap_Hz = get_objectBasedFieldmap_Hz(signal, voxelSize_mm, B0dir, centerFreq_Hz);
signal = demodulate_field_Hz(signal, TE_s, objectBasedFieldmap_Hz);
p3 = angle(signal(:, :, isl, end));

% residual linear field
cropDim = 1;
residualLinearField_Hz = get_residualLinearField_Hz(signal, TE_s, cropDim);

signal = demodulate_field_Hz(signal, TE_s, residualLinearField_Hz);
p4 = angle(signal(:, :, isl, end));



%% water--fat separation

imDataParams.images = reshape(signal(:, :, isl, :), ...
                              [sz(1), sz(2), 1, 1, sz(4)]); % add dummy coil dimension, TE dimension last
imDataParams.TE = TE_s;
imDataParams.FieldStrength = 3;
imDataParams.PrecessionIsClockwise = 1; % Philips scanner follows convention of clockwise precession

algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat (7 peaks)';
algoParams.species(2).frequency = [-3.8, -3.4, -3.1, -2.68, -2.46, -1.95, -0.5, 0.49, 0.59];
algoParams.species(2).relAmps = [0.0899, 0.5834, 0.0599, 0.0849, 0.0599, 0.0150, 0.0400, 0.01, 0.0569];


% hIDEAL
algoParams.Visualize = 0;
algoParams.AlwaysShowGUI = 0;

outParams_hIDEAL = fw_i2cm0c_3pluspoint_tsaojiang(imDataParams, algoParams);

% graph cut
algoParams.range_fm = [-1000, 1000];% Range of field map values
algoParams.NUM_FMS = 500;% Number of field map values to discretize

outParams_graphcut = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams, algoParams);


% without demodulation steps for comparison
imDataParams.images = reshape(signal0(:, :, isl, :), ...
                              [sz(1), sz(2), 1, 1, sz(4)]); % add dummy coil dimension, TE dimension last
outParams0_hIDEAL = fw_i2cm0c_3pluspoint_tsaojiang(imDataParams, algoParams);
outParams0_graphcut = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams, algoParams);



%% plot results

% standard WFI results withouth any demodulation steps
figure('position', [0, 0, 1000, 1000])
colormap gray
subplot(2, 2, 1)
imagesc(abs(outParams0_hIDEAL.species(1).amps))
set(gca, 'xtick', [])
set(gca, 'ytick', [])
ylabel('water')
title('hIDEAL')
subplot(2, 2, 2)
imagesc(abs(outParams0_graphcut.species(1).amps))
axis off
title('graph cut')
subplot(2, 2, 3)
imagesc(abs(outParams0_hIDEAL.species(2).amps))
set(gca, 'xtick', [])
set(gca, 'ytick', [])
ylabel('fat')
subplot(2, 2, 4)
imagesc(abs(outParams0_graphcut.species(2).amps))
axis off

% WFI results after proposed demodulation steps
figure('position', [1000, 0, 1000, 1000])
colormap gray
subplot(2, 2, 1)
imagesc(abs(outParams_hIDEAL.species(1).amps))
set(gca, 'xtick', [])
set(gca, 'ytick', [])
ylabel('water')
title('hIDEAL')
subplot(2, 2, 2)
imagesc(abs(outParams_graphcut.species(1).amps))
axis off
title('graph cut')
subplot(2, 2, 3)
imagesc(abs(outParams_hIDEAL.species(2).amps))
set(gca, 'xtick', [])
set(gca, 'ytick', [])
ylabel('fat')
subplot(2, 2, 4)
imagesc(abs(outParams_graphcut.species(2).amps))
axis off

% row 1: magnitude image TE1, estimated field contributions in [Hz]
% row 2: phase images of TE3 after demodulation of field contribution in row 1
figure('position', [0, 0, 2500, 1000])
colormap('gray')
subplot(2, 5, 1)
imagesc(abs(signal(:, :, isl, 1)))
axis off
title('magnitude TE1')
subplot(2, 5, 2)
imagesc(magnetInhomogeneities_Hz(:, :, isl))
colorbar
axis off
title('B0 inhom.')
subplot(2, 5, 3)
imagesc(shimField_Hz(:, :, isl))
title('shim')
colorbar
axis off
subplot(2, 5, 4)
imagesc(objectBasedFieldmap_Hz(:, :, isl))
title('OBFFME')
colorbar
axis off
subplot(2, 5, 5)
imagesc(residualLinearField_Hz(:, :, isl))
title('res. lin. field.')
colorbar
axis off

subplot(2, 5, 6)
imagesc(p0)
ylabel('phase TE3')
set(gca, 'xtick', [])
set(gca, 'ytick', [])
subplot(2, 5, 7)
imagesc(p1)
axis off
subplot(2, 5, 8)
imagesc(p2)
axis off
subplot(2, 5, 9)
imagesc(p3)
axis off
subplot(2, 5, 10)
imagesc(p4)
axis off
