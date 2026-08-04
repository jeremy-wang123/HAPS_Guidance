function results = estimate_insar_uncertainty( ...
    L, lambda, dphi, NESZ_dB, sigma0_dB, b_perp, time)
%ESTIMATE_INSAR_UNCERTAINTY Estimate InSAR coherence and uncertainty.
%
% Inputs
%   L           - Number of independent looks
%   lambda      - Radar wavelength, m
%   dphi        - Difference in platform heading, degrees
%   NESZ_dB     - Noise Equivalent Sigma Zero, dB
%   sigma0_dB   - backscatter coefficient, dB
%   b_perp      - Perpendicular baseline, m
%   time        - Temporal separation, seconds; may be scalar or array
%
% Outputs
%   results - Structure containing:
%       .time
%       .gamma_temporal
%       .gamma_thermal
%       .gamma_spatial
%       .gamma_rotation
%       .gamma_total
%       .sigma_phi
%       .sigma_z
%
% Fixed system parameters in this function
%   theta   = 40 degrees
%   ry      = ground-range resolution
%   rx      = azimuth resolution
%   altitude = 18 km
%   T_zero   = 6 days
%
% Notes
%   1. Temporal coherence decreases linearly to zero at six days.
%   2. Coherence values are constrained to the interval [0,1].
%   3. sigma_phi and sigma_z are undefined when total coherence is zero.

%% Input validation

validateattributes(L, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'L');

validateattributes(lambda, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'}, mfilename, 'lambda');

validateattributes(dphi, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, 'dphi');

validateattributes(NESZ_dB, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, 'NESZ_dB');

validateattributes(sigma0_dB, {'numeric'}, ...
    {'scalar', 'real', 'finite'}, mfilename, 'sigma0_dB');

validateattributes(b_perp, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, 'b_perp');

validateattributes(time, {'numeric'}, ...
    {'real', 'finite', 'nonnegative'}, mfilename, 'time');

%% Fixed radar and acquisition parameters

theta = 40;                 % look angle, degrees
r_slant = 0.75;             % slant-range resolution, m
rx = 0.5;                   % azimuth resolution, m
altitude = 18e3;            % platform altitude, m

% Ground-range resolution
ry = r_slant / sind(theta);

% Slant range from platform to ground
R = altitude / cosd(theta);

% Time at which temporal coherence reaches zero
T_zero = 6 * 24 * 3600;     % seconds

%% Temporal coherence

% Linear decrease from one at time = 0 to zero at time = T_zero
gamma_temporal = max(1 - time ./ T_zero, 0);

%% Thermal coherence
SNR_dB = sigma0_dB - NESZ_dB;
SNR = 10^(SNR_dB/10);

% SNR must be in linear units
gamma_thermal = SNR / (SNR + 1);

%% Spatial coherence

% Critical perpendicular baseline
Bc = (lambda * R) / (2 * ry * cosd(theta));

if b_perp >= Bc
    gamma_spatial = 0; 
else
    gamma_spatial = 1 - (2 * abs(b_perp) * ry * cosd(theta)) / (lambda * R);
end


%% Rotational coherence

dphi_rad = deg2rad(dphi);

dphi_c = rad2deg(lambda / (2*rx * sind(theta))); % maximum heading difference allowed

if abs(dphi) >= dphi_c
    gamma_rotation = 0;
else
    gamma_rotation = 1 - (2*rx*dphi_rad*sind(theta))/lambda; % rotational coherence
end


%% Total coherence

gamma_total = ...
    gamma_thermal .* ...
    gamma_spatial .* ...
    gamma_rotation .* ...
    gamma_temporal;

gamma_total = max(min(gamma_total, 1), 0);

%% Phase uncertainty

sigma_phi = NaN(size(gamma_total));

% Phase uncertainty is only finite when total coherence is greater than zero
valid = gamma_total > 0;

sigma_phi(valid) = ...
    sqrt(1 - gamma_total(valid).^2) ./ ...
    (gamma_total(valid) .* sqrt(2 * L));

%% Vertical displacement uncertainty

% Assumes displacement is entirely vertical:
%
% delta_r_LOS = delta_z*cos(theta)
%
% sigma_z = lambda/(4*pi*cos(theta))*sigma_phi

sigma_z = NaN(size(sigma_phi));

sigma_z(valid) = ...
    (lambda / (4 * pi * cosd(theta))) .* sigma_phi(valid);

%% Store outputs

results.time = time;
results.gamma_temporal = gamma_temporal;

results.gamma_thermal = gamma_thermal;
results.gamma_spatial = gamma_spatial;
results.gamma_rotation = gamma_rotation;
results.gamma_total = gamma_total;

results.sigma_phi = sigma_phi;
results.sigma_z = sigma_z;

% Useful intermediate quantities
results.theta = theta;
results.rx = rx;
results.ry = ry;
results.R = R;
results.Bc = Bc;
results.dphi_c = dphi_c;
results.T_zero = T_zero;

end