function img = shot_gather_mig(x, z, c_x_z, f, rxdata_f, txdata_f, aawin)

% Verify the Number of Common Shot Gathers
ns = size(txdata_f, 3); assert(size(rxdata_f, 3) == ns, ...
    'Number of sources must equal to number of common-source gathers');
AAwin = repmat(aawin(:), [1, ns]);

% Forward and Inverse Fourier Transforms with Anti-Aliasing Windows
ft = @(sig) fftshift(fft(AAwin.*sig, [], 1), 1);
ift = @(sig) AAwin.*ifft(ifftshift(sig, 1), [], 1);

% Spatial Grid
dx = mean(diff(x)); nx = numel(x); 
x = dx*((-(nx-1)/2):((nx-1)/2)); 
dz = mean(diff(z)); 

% FFT Axis for Lateral Spatial Frequency
kx = mod(fftshift((0:nx-1)/(dx*nx))+1/(2*dx), 1/dx)-1/(2*dx);

% Convert to Slowness [s/m]
s_x_z = 1./c_x_z; % Slowness = 1/(Speed of Sound)
s_z = mean(s_x_z, 2); % Mean Slowness vs Depth (z)
ds_x_z = s_x_z - repmat(s_z, [1, numel(x)]); % Slowness Deviation

% Construct Image by Looping Over Frequency
img = zeros(numel(z), numel(x));
M = moviein(numel(f));
for f_idx = 1:numel(f)
    % Setup Downward Continuation at this Frequency
    rx_surface = squeeze(permute(rxdata_f(:,f_idx,:), [2,1,3])); 
    tx_surface = squeeze(permute(txdata_f(:,f_idx,:), [2,1,3])); 
    rx_singleFreq_x_z = zeros(numel(z), numel(x), ns);
    tx_singleFreq_x_z = zeros(numel(z), numel(x), ns);
    rx_singleFreq_x_z(1, :, :) = rx_surface;
    tx_singleFreq_x_z(1, :, :) = tx_surface;
    % Continuous Wave Response By Downward Angular Spectrum
    for z_idx = 1:numel(z)-1
        % Create Propagation Filter for this Depth
        kz = sqrt((f(f_idx)*s_z(z_idx)).^2 - kx.^2); % Axial Spatial Frequency
        H = exp(1i*2*pi*kz*dz); % Propagation Filter in Spatial Frequency Domain
        H(kz.^2 <= 0) = 0; % Remove Evanescent Components
        H = repmat(H(:), [1, ns]); % Replicate Across Shots
        % Create Phase-Shift Correction in Spatial Domain
        dH = exp(1i*2*pi*f(f_idx)*ds_x_z(z_idx,:)*dz); 
        dH = repmat(dH(:), [1, ns]); % Replicate Across Shots
        % Downward Continuation with Split-Stepping
        rx_singleFreq_x_z(z_idx+1, :, :) = dH .* ...
            ift(H.*ft(reshape(rx_singleFreq_x_z(z_idx,:,:),[numel(x),ns])));
        tx_singleFreq_x_z(z_idx+1, :, :) = conj(dH) .* ...
            ift(conj(H).*ft(reshape(tx_singleFreq_x_z(z_idx,:,:),[numel(x),ns]))); 
    end
    % Accumulate Image Frequency-by-Frequency
    img = img + sum(tx_singleFreq_x_z .* conj(rx_singleFreq_x_z), 3);
    % Reconstruct and Plot Ultrasound Image
    imagesc(1000*x, 1000*z, db(abs(img)/max(abs(img(:)))), [-80, 0]);
    xlabel('x Azimuthal Distance (mm)'); ylabel('z Axial Distance (mm)'); 
    title(['Image Reconstruction up to ', num2str(f(f_idx)/(1e6)), ' MHz']); 
    zoom on; axis equal; axis xy; axis image; colormap gray; colorbar(); 
    xlim([-19.1,19.1]); set(gca, 'YDir', 'reverse'); M(f_idx) = getframe(gcf);
end
% Save Accumulation of Image in Frequency Domain
%movie2gif(M, 'FreqDomain.gif');

end

