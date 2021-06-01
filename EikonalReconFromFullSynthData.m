clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));

% Load Channel Data and Sound Speed Estimate
load ChickenPhantomMultiFocal.mat; % Channel Data
load SoundSpeedEstimate.mat; % Estimated Sound Speed Map

% Select Whether or Not You Want to Correct for Aberration
aberration_corr = 1; % 1 for aberration correction, 0 for regular DAS

% Select Subset of Transmit Elements
no_elements = size(rxAptPos,1);
tx_elmts = 1:4:no_elements;
txAptPos = rxAptPos(tx_elmts,:);
rxdata = scat(:,:,tx_elmts);
rxdata_h = reshape(hilbert(reshape(rxdata, ...
    [numel(time), no_elements*numel(tx_elmts)])), ...
    [numel(time), no_elements, numel(tx_elmts)]);

% Points to Focus and Get Image At
num_x = 301; xlims = (20e-3)*[-1, 1];
num_z = 601; zlims = [0e-3, 45e-3];
x_img = linspace(xlims(1), xlims(2), num_x);
z_img = linspace(zlims(1), zlims(2), num_z);
dBrange = [-80, 0]; c_liver = 1540;

% Arrival Times
[X, Z] = meshgrid(x, z); dx = mean(diff(x));
[X_img, Z_img] = meshgrid(x_img, z_img);
foc_pts = [X_img(:), 0*Z_img(:), Z_img(:)];
t_rx = zeros(num_z, num_x, no_elements);
t_rx_c = zeros(num_z, num_x, no_elements);
for elmt = 1:no_elements
    % Calculating Source Point Location
    [~, Iz] = min(abs(z-rxAptPos(elmt,3)));
    [~, Ix] = min(abs(x-rxAptPos(elmt,1)));
    % Travel Time Calculation
    t_tx = dx*msfm2d(Crecon, [Iz; Ix], true, true); disp(elmt);
    % Interpolation onto Imaging Grid
    t_rx(:,:,elmt) = interp2(X, Z, t_tx, X_img, Z_img, 'spline');
    t_rx_c(:,:,elmt) = reshape(calc_times(foc_pts,...
        rxAptPos(elmt,:),c_liver),[num_z,num_x]);
end

% Full Synthetic Aperture Image Reconstruction
if aberration_corr
    rx_times = reshape(t_rx, [num_z*num_x, no_elements]); 
else
    rx_times = reshape(t_rx_c, [num_z*num_x, no_elements]); 
end
tx_times = rx_times(:,tx_elmts);
tic; focData = focus_eikonal(time, ...
    rxdata_h, rx_times, tx_times); toc;
focData = reshape(focData, [num_z, num_x]); 
img_h = sum(sum(focData,3),4);

% Two Methods for Displaying the Ultrasound Images
imagesc(1000*x_img, 1000*z_img, ...
    20*log10(abs(img_h)/max(abs(img_h(:)))), dBrange); 
axis image; xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title('Eikonal DAS Beamforming'); colormap(gray); colorbar();