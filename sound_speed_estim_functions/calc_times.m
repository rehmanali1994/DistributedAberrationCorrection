function foc_times = calc_times(foci, elempos, speed_of_sound)
% foc_times = calc_times(foci, elempos, speed_of_sound)
% 
% calc_times - computes focusing times
% 
% The function computes the (Tx or Rx) time of arrival for specified focal points
% given the array element positions.
% 
% NOTE: Primarily intended when Tx and Rx apertures are the same (i.e. no full synthetic aperture)
% 
% INPUTS:
% foci              - M x 3 matrix with position of focal/image points of interest [m]
% elempos           - N x 3 matrix with element positions [m]
% speed_of_sound    - speed of sounds [m/s]
% 
% OUTPUT:
% foc_times         - M x N matrix with travel times between each image point and array element

foci_tmp = repmat(reshape(foci, ...
    [size(foci,1),1,size(foci,2)]),[1,size(elempos,1),1]);
elempos_tmp = repmat(reshape(elempos, ...
    [1,size(elempos,1),size(elempos,2)]), [size(foci,1),1,1]);
distance = sqrt(sum((foci_tmp-elempos_tmp).^2,3));
foc_times = distance/speed_of_sound;

end