function A = linRayTimesSparseMatrix(x, z, soundSpeedMap, tx_loc, rx_loc)
%LINRAYTIMESSPARSE Calculate Arrival Times in Space from Each Array
%Element From Sound Speed Map Using a Straight Ray Path Approximation 
%   x -- Lateral Coordinate [m] as 1 x Nx vector
%   z -- Axial Coordinate[m] as 1 x Nz vector
%   soundSpeedMap -- Sound Speed Map [m/s] as M x N Array
%   tx_loc -- (x, z) location of transmit elements as 2 x Ktx Array
%   rx_loc -- (x, z) location of receive elements as 2 x Krx Array
%   A -- Sparse Tomography Matrix (Ktx*Krx x Nz*Nx Array) 

% Calculate Slowness
slownessMap = 1./soundSpeedMap;

% Grid Sampling
dx = mean(diff(x)); dz = mean(diff(z));
Nx = numel(x); Nz = numel(z);
[X, Z] = meshgrid(x, z);

% Ray Sampling
dr = 0.1/((1/dx)+(1/dz));

% Create Arrival Time Array for Polar Coordinates
Ktx = size(tx_loc, 2); Krx = size(rx_loc, 2);
A = sparse(Ktx*Krx, Nz*Nx);

% Build Sparse System Matrix for Straight-Path Forward Projection
for tx_idx = 1:Ktx
    % Calculate Eikonal Solution and Ray Tracing
    T = sqrt((X-tx_loc(1,tx_idx)).^2+(Z-tx_loc(2,tx_idx)).^2);
    tic; [path_x, path_z] = rayTracEikonDescnd(x, z, T, rx_loc(1,:), ...
        rx_loc(2,:), tx_loc(1,tx_idx), tx_loc(2,tx_idx), dr); toc;
    % Accumulate Sparse Matrix for these Projections
    tic; path_x_idx = round(min(Nx, max(1, 1+(path_x-x(1))/dx)));
    path_z_idx = round(min(Nz, max(1, 1+(path_z-z(1))/dz)));
    path_x_idx(isnan(path_x)) = NaN;
    path_z_idx(isnan(path_z)) = NaN;
    path_ind = sub2ind([Nz, Nx], path_z_idx, path_x_idx);
    elmt_ind = ones(size(path_ind,1),1)*(1:Krx)+(tx_idx-1)*Krx;
    elmt_ind_flat = elmt_ind(~isnan(path_ind));
    path_ind_flat = path_ind(~isnan(path_ind));
    A = A + dr*sparse(elmt_ind_flat, path_ind_flat, ...
        ones(numel(elmt_ind_flat),1), Ktx*Krx, Nz*Nx); toc;
    disp(['Sparse Matrix Tx Element ', num2str(tx_idx), ' Completed']);
end

end