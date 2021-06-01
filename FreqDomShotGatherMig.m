clear
clc

% Load all Functions from Subdirectories
addpath(genpath(pwd));

% Load Channel Data and Sound Speed Estimate
load ChickenPhantomMultiFocal.mat; % Channel Data
load SoundSpeedEstimate.mat; % Estimated Sound Speed Map

% Set Imaging Grid
dov = 45e-3; % Max Depth [m]
upsamp_x = 2; % Upsampling in x (assuming dx = pitch)
upsamp_z = 2; % Upsampling in z (assuming dz = pitch)
Nx0 = 256; % Number of Points Laterally in x
[X,Z] = meshgrid(x,z); % Coordinate for Reconstructed Sound Speed

% Pick off Specific Transmitter
[nT, nRx, nTx] = size(scat); 
scat_h = reshape(hilbert(reshape(scat, ...
    [nT,nRx*nTx])),[nT,nRx,nTx]);
tx_elmts = 1:2:no_elements; clearvars scat; 
rxdata_h = scat_h(:,:,tx_elmts);

% Aperture Definition
xpos = rxAptPos(:,1); 
pitch = mean(diff(xpos)); 
no_elements = numel(xpos);

% Simulation Space
x = (-(upsamp_x*Nx0-1)/2:(upsamp_x*Nx0-1)/2)*(pitch/upsamp_x); % m
Nu1 = round(dov/(pitch/upsamp_z)); 
z = ((0:Nu1-1))*pitch/upsamp_z;
[Xg, Zg] = meshgrid(x,z); 
c_x_z = interp2(X, Z, Crecon, Xg, Zg);

% Some Logic to Make Sound Speed Continuous Outside Estimated Domain
interior_idx = find(~any(isnan(c_x_z)));
c_left = c_x_z(:, interior_idx(1));
c_right = c_x_z(:, interior_idx(end));
c_x_z(:, any(isnan(c_x_z)) & (x < 0)) = ...
    repmat(c_left, [1, sum(any(isnan(c_x_z)) & (x < 0))]); 
c_x_z(:, any(isnan(c_x_z)) & (x > 0)) = ...
    repmat(c_right, [1, sum(any(isnan(c_x_z)) & (x > 0))]); 

% Image Reconstruction Parameters and Anti-Aliasing Window
dBrange = [-80, 0]; ord = 50; 
xmax = (max(abs(xpos))+max(abs(x)))/2;
aawin = 1./sqrt(1+(x/xmax).^ord);

% Time and Frequency Axes
nt = numel(time); fs = 1/mean(diff(time));
f = (fs/2)*(-1:2/nt:1-2/nt); % Hz

% Get Receive Channel Data in the Frequency Domain
P_Rx_f = fftshift(fft(rxdata_h, nt, 1), 1);
[~, F, ~] = meshgrid(1:size(rxdata_h,2), f, 1:size(rxdata_h,3)); 
P_Rx_f = P_Rx_f .* exp(-1i*2*pi*F*time(1));
rxdata_f = interp1(xpos, permute(P_Rx_f, [2,1,3]), x, 'nearest', 0);

% Only Keep Positive Frequencies within Transducer Passband
passband_f_idx = find((f > 3e6) & (f < 12e6));
rxdata_f = rxdata_f(:,passband_f_idx,:);
P_Tx_f = ones(size(f)); f = f(passband_f_idx); 
P_Tx_f = P_Tx_f(passband_f_idx); % Assume Flat Passband

% Pulsed-Wave Frequency Response on Transmit
apod = eye(no_elements); delay = zeros(no_elements);
txdata_f = zeros(numel(x), numel(f), numel(tx_elmts));
for k = 1:numel(tx_elmts) 
    % Construct Transmit Responses for Each Element
    apod_x = interp1(xpos, apod(:,tx_elmts(k)), x, 'nearest', 0);
    delayIdeal = interp1(xpos, delay(:,tx_elmts(k)), x, 'nearest', 0);
    txdata_f(:,:,k) = (apod_x'*P_Tx_f).*exp(-1i*2*pi*delayIdeal'*f);
end

% Reconstruct Ultrasound Image Using Shot-Gather Migration
img = shot_gather_mig(x, z, c_x_z, f, rxdata_f, txdata_f, aawin);

% Generate Time-Gain Compensation From Center Frequency
f_ctr_idx = round(numel(f)/2);
tx_map = zeros(numel(z), numel(x));
for elmt = 1:numel(tx_elmts)
    tx_singleFreq_x_z = angular_spectrum(x, z, c_x_z, ...
        f(f_ctr_idx), txdata_f(:,f_ctr_idx,elmt), aawin);
    tx_map = tx_map + tx_singleFreq_x_z .* conj(tx_singleFreq_x_z);
end

% Apply Time-Gain Compensation to Image
reg = 1e-2; img_recon = img ./ (tx_map + reg*max(tx_map(:)));

% Reconstruct and Plot Ultrasound Image
imagesc(1000*x, 1000*z, db(abs(img_recon)/max(abs(img_recon(:)))), dBrange);
xlabel('x Azimuthal Distance (mm)'); ylabel('z Axial Distance (mm)'); 
title('Image Reconstruction'); zoom on; axis equal; axis xy; axis image; 
colormap gray; colorbar(); set(gca, 'YDir', 'reverse'); xlim([-19.1,19.1]);
