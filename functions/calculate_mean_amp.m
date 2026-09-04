function mean_amp = calculate_mean_amp(glacier_name, top_n)

% CALCULATE_MEAN_AMP calculates the mean summed tidal amplitude
% over the floating ice shelf for the top N tidal constituents.
%
% Inputs:
%   glacier_name - Glacier name, e.g. 'Thwaites'
%   top_n        - Number of top constituents to use (1-10)
%
% Output:
%   mean_amp     - Mean summed tidal amplitude over floating ice shelf
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008')

Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008';

boundary_name = glacier_name;
region_name   = [glacier_name ' Glacier'];

%% Check input
if top_n < 1 || top_n > 10 || mod(top_n,1) ~= 0
    error('top_n must be an integer between 1 and 10.');
end

%% Extract grid
[xgrid,ygrid] = psgrid(region_name,700,1,'xy');

%% Get constituent amplitudes
% Get all constituents ranked by mean amplitude
[constit_list, constit_mean_amp, constit_idx] = ...
    select_constituent_amplitudes(region_name, Model, 'all');

% Keep only top N
constit_idx = constit_idx(1:top_n);

%% Extract tidal amplitudes
[lat,lon] = ps2ll(xgrid,ygrid);

[amps,~,~,~] = tmd_extract_HC( ...
    Model, lat, lon, 'z', constit_idx);

if top_n == 1
    amp_sum = amps;
else
    % Sum across tidal constituents
    amp_sum = squeeze(sum(amps,1));
end
%% Get glacier boundary
[wx,wy] = antbounds_data(boundary_name,'xy');

inside_region = inpolygon(xgrid,ygrid,wx,wy);

%% Identify floating ice
shelf = isiceshelf(xgrid,ygrid);

glacier_shelf = inside_region & shelf;

%% Calculate mean amplitude
mean_amp = mean(amp_sum(glacier_shelf),'omitnan');

end