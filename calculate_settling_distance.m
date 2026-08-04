function [settling_distance] = calculate_settling_distance(w, x, w0, sigma_z)
%CALCULATE_SETTLING_DISTANCE calculates the distance at which the beam
%profile settles within the range of displacement uncertainty
% Inputs:
%   - w (w(x) profile)
%   - x (x space)
%   - w0 (tidal amplitude)
%   - sigma_z (vertical displacement uncertainty)
% Outputs:
%   - settling_distance (km)
    residual = w - w0;
    
    outside_band = abs(residual) > sigma_z;
    
    last_outside = find(outside_band, 1, 'last');
    
    if isempty(last_outside)
    
        % Entire profile is already within the uncertainty band
        settling_distance_m = 0;
    
    elseif last_outside == numel(w)
    
        % Domain is not long enough to observe settling
        settling_distance_m = NaN;
    
        warning(['The profile has not remained within +/- sigma_z ', ...
                 'by the end of the modeled domain.']);
    
    else
    
        % First point after the final departure from the uncertainty band
        settling_distance_m = x(last_outside + 1);
    
    end
    
    settling_distance = settling_distance_m/1000; % km
end
