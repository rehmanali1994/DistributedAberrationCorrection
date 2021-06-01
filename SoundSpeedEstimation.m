clear
clc

% Load Saved CF Images
load('SavedCoherenceFactorImages.mat');
load('ChickenPhantomMultiFocal.mat', 'rxAptPos');

%%% 1) Average Sound Speed Estimates

% Assemble Coherence Factor Curves for Each Depth
u = gausswin(100, 2); % Filtering in Axial Dimension
v = gausswin(50, 2); % Filtering in Lateral Dimension
cf_imgs = zeros(size(cf)); % Blurred CF Images

% Filter All Metrics
for k = 1:numel(c_bfm)
    cf_imgs(:,:,k) = conv2(u, v, reshape(cf(:,:,k), ...
        [numel(z_img), numel(x_img)]), 'same');
end

% Estimate Average Sound Speed up to Imaging Depth
c_bfm_increment = 0.1; % Upsampling Increment for Sound Speed Estimation
c_bfm_upsamp = c_bfm(1):c_bfm_increment:c_bfm(end);
cf_curves_upsamp = interp1(c_bfm, permute(cf_imgs,[3,1,2]), c_bfm_upsamp);
[~,index] = max(cf_curves_upsamp); 
c_avg = squeeze(c_bfm_upsamp(index));

% Smoothing Regularization
reg = 1.0e7; % Regularization Value
I = eye(numel(z_img)); % Identity Matrix
F = I(1:end-3,:)-3*I(2:end-2,:)+3*I(3:end-1,:)-I(4:end,:);
c_avg_smoothed = (I+reg*F'*F)\c_avg;

%%% 2) Setup for Local Sound Speed Estimation

% Create Grid for Sound Speed Map
z_grid_pad = 10; % Extend Grid Behind Transducer
no_subelem = 3; % Number of Mathematical Subelements per Element
pitch = mean(diff(rxAptPos(:,1))); % Pitch [m]
dx = pitch/no_subelem; % Grid Spacing [m]
x = (rxAptPos(1,1)-dx*(no_subelem-1)/2):dx:...
    (rxAptPos(end,1)+dx*(no_subelem-1)/2); % Lateral Spatial Grid
z = -z_grid_pad*dx:dx:z_img(end); % Axial Spatial Grid
Nx = numel(x); Nz = numel(z); % Number of Grid Points in x and z
[X, Z] = meshgrid(x, z); % Create Meshgrid
c_guess = 1540; % Initial Sound Speed [m/s] Guess 
C = c_guess*ones(size(X)); % Sound Speed [m/s] Map

% Receiver Locations
x_idx = 5:5:145; xrx = x_img(x_idx); Nx_rx = numel(xrx);
z_idx = 5:10:595; zrx = z_img(z_idx); Nz_rx = numel(zrx);
[Xrx, Zrx] = meshgrid(xrx, zrx);
xrx = Xrx(:)'; zrx = Zrx(:)';
rx_loc = [xrx; zrx];

% Transmitter Locations
elmt_skip = 2;
xtx = rxAptPos(1:elmt_skip:end, 1)';
ztx = rxAptPos(1:elmt_skip:end, 3)';
Nelem = numel(xtx);
tx_loc = [xtx; ztx]; 

% Create Tomography Matrix for Sound Speed Reconstruction
if ~exist('tomography_matrix.mat','file')
    A = linRayTimesSparseMatrix(x, z, C, tx_loc, rx_loc);
    save('tomography_matrix.mat', ...
        '-v7.3','x','z','C','tx_loc','rx_loc','A');
end
load('tomography_matrix.mat','A'); D = A; clearvars A;

% Convert Average Sound Speed Measurements into Travel Times
c_avg = c_avg_smoothed(z_idx, x_idx); % Subsampled Image Grid
times = zeros(Nz_rx, Nx_rx, Nelem);
for x_idx = 1:Nx_rx
    for z_idx = 1:Nz_rx
        % Convert Sound Speed Fit into Travel Times
        distance = sqrt((xtx-Xrx(z_idx,x_idx)).^2+(ztx-Zrx(z_idx,x_idx)).^2);
        times(z_idx,x_idx,:) = distance/c_avg(z_idx,x_idx);
    end
end

%%% 3) Local Sound Speed Estimation

% Dimensions of H
M = Nz_rx*Nx_rx*Nelem; 
N = Nz*Nx;

% Weighting
[Xtx, Xrx] = meshgrid(xtx, xrx); 
[Ztx, Zrx] = meshgrid(zeros(Nelem,1), zrx); 
Dsq = (Xtx(:)-Xrx(:)).^2 + (Ztx(:)-Zrx(:)).^2; 
W = (Dsq(:)/max(Dsq(:))).^(1);

% Forward Model Matrix as Linear Operator
H = @(x) D*x;
HT = @(x) D'*x;

% Prior Covariance Matrix Q as a Linear Operator
u = gausswin(200, 4); v = gausswin(500, 4);
Q = @(x) reshape(conv2(u, v, reshape(x, [Nz, Nx]), 'same'), [Nz*Nx, 1]);

% Observation Covariance Matrix R as Sparse Matrix
R = (1e-2)*W;

% Prior Mean Vector
c_guess = 1540; % Uniform Sound Speed Guess
mu = (1/c_guess)*ones(size(C(:))); % Running Reconstruction

% Define A (Symmetric) Matrix as Linear Operator
A = @(x) H(Q(HT(x)))+R.*x;

% Calculate Vector b
y = times(:);
b = y - H(mu);

% Initial Conditions for Conjugate Gradient Algorithm
xi = zeros(M, 1); % Initial Solution
beta = 0; % Variable for Updating Conjugate Direction
p = zeros(size(xi)); % Conjugate Direction Vector
r = A(xi)-b; % Current Gradient Direction

% Iterations of CG
num_iter = 400; nskip = 10;
% Show Average Sound Speed Measurements
subplot(1,2,1); imagesc(1000*x_img, 1000*z_img, c_avg_smoothed); 
axis image; colorbar; hold on; caxis([1450,1650])
xlim(1000*[min(x), max(x)]); ylim(1000*[min(z), max(z)]); 
plot(1000*xtx, 1000*ztx, 'k.', 'Linewidth', 1); 
plot(1000*xrx, 1000*zrx, 'w.', 'Linewidth', 1); 
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); 
title('Estim. Avg Sound Speed [m/s]');
% Iterations
for k = 1:num_iter
    % Conjugate Gradient Updates
    p = -r + beta*p;
    r_norm_sq_last = r'*r; Ap = A(p);
    alpha = r_norm_sq_last/(p'*Ap);
    xi = xi + alpha*p;
    r = r + alpha*Ap;
    beta = (r'*r)/r_norm_sq_last;  
    if mod(k, nskip) == 1
        % Calculate Posterior Mean (Reconstructed Image
        mupost = mu + Q(HT(xi));
        % Show Reconstructed Image of Sound Speed Map
        subplot(1,2,2); Crecon = reshape(1./mupost, [Nz, Nx]);
        imagesc(1000*x, 1000*z, Crecon, [1450, 1650]); 
        colorbar; axis image; hold on;
        plot(1000*xtx, 1000*ztx, 'k.', 'Linewidth', 1);
        plot(1000*xrx, 1000*zrx, 'w.', 'Linewidth', 1);
        xlabel('Axial [mm]'); ylabel('Lateral [mm]');
        title(['Conjugate Gradient Iteration ', num2str(k)]);
        getframe(gca);
    end
end

% Show Image of Final Reconstructed Sound Speed Map
subplot(1,2,2); Crecon = reshape(1./mupost, [Nz, Nx]);
imagesc(1000*x, 1000*z, Crecon, [1450, 1650]); 
colorbar; axis image; hold on; 
plot(1000*xtx, 1000*ztx, 'k.', 'Linewidth', 1);
plot(1000*xrx, 1000*zrx, 'w.', 'Linewidth', 1);
ylabel('Axial [mm]'); xlabel('Lateral [mm]');
title('Tomographic Sound Speed Recon. [m/s]');

% Compute Travel Times from Slowness
slowness = 1./C; t = D*mupost;
t = reshape(t, [Nz_rx, Nx_rx, Nelem]);

% Create Grid on Which to Calculate Travel Times via Eikonal Equation
Crecon = Crecon(z_grid_pad+1:end,:);
z = z(z_grid_pad+1:end); 

% Save Estimated Sound Speed Reconstruction
save('SoundSpeedEstimate','-v7.3','x','z','Crecon');