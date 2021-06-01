function cf_img = multistatic_cf(t, signal, foc_pts, rxAptPos, txAptPos, speed_of_sound)
% cf_img = multistatic_cf(t, signal, foc_pts, rxAptPos, txAptPos, speed_of_sound)
% 
% multistatic_cf - Coherence Factor (CF) Beamforming of Multistatic Focused RF Data
% 
% The function interpolates the RF signals collected using the full synthetic sequence
% to focus the data at desired locations.
% 
% INPUTS:
% t                  - T x 1 time vector for samples of the input signal [s]
% signal             - T x N x M matrix containing input RF data to be interpolated
% foc_pts            - P x 3 matrix with position of focal points [m]
% rxAptPos           - N x 3 matrix with positions of the Rx elements [m]
% txAptPos           - M x 3 matrix with positions of the Tx elements [m]
% speed_of_sound     - speed of sounds [m/s]
% 
% OUTPUT:
% cf_img             - vector with dimension P for CF image

% time from the focus to receive  apertures (array elements)
rx_times = calc_times(foc_pts, rxAptPos, speed_of_sound);

% time from the transmit apertures (array elements) to focus
tx_times = calc_times(foc_pts, txAptPos, speed_of_sound);

% compute coherent and incoherent focused images
coherent_img = zeros(size(foc_pts,1),1);
incoherent_img = zeros(size(foc_pts,1),1);
for i = 1:size(rx_times,2)
    % complex image for single receive element
    rx_img = zeros(size(foc_pts,1),1);
    for j = 1:size(tx_times,2)
        rx_img = rx_img + interp1(t, signal(:,i,j), ...
            rx_times(:,i)+tx_times(:,j), 'spline', 0);
    end
    % Compute Coherent and Incoherent Sums
    coherent_img = coherent_img + rx_img; 
    incoherent_img = incoherent_img + abs(rx_img).^2;
end

% compute coherence factor images
cf_img = (abs(coherent_img).^2)./(incoherent_img*size(rx_times,2));
    
end