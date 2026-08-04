function [h_shelf_mean,h_shelf_eff,h_front_mean,h_front_eff] = extract_shelf_thickness(region)
%EXTRACT_SHELF_THICKNESS extracts the mean thickness and calving front thickness of a
%given ice shelf
%   Calving front thickness averaged across data points within 5 km of the
%   calving front. 
%   Input: region (i.e. Thwaites Glacier)
%   Outputs: 
%   - h_shelf_mean: average thickness of ice shelf
%   - h_shelf_eff: rigidity equivalent mean thickness
%   - h_front_mean: average thickness of region within 5km of calving front
%   - h_front_eff: rigidity equivalent thickness of calving front
%   Data set from Bedmachine (Morlighem et al. 2017)

addpath('/Users/jeremywang/Documents/MATLAB/BedMachine');

front_band_width_km = 5; % Near-front region used for summary statistics [km]
grid_spacing_km = 1;
grid_spacing_m = grid_spacing_km*1000;

% extract coordinates
[xgrid,ygrid] = psgrid('Thwaites Glacier',700,grid_spacing_km,'xy');

% Extract BedMachine data
H = bedmachine_interp( ...
    'thickness', ...
    xgrid, ...
    ygrid);

bm_mask = bedmachine_interp( ...
    'mask', ...
    xgrid, ...
    ygrid, ...
    'method', ...
    'nearest');

% BedMachine mask values:
%   0 = ocean
%   1 = ice-free land
%   2 = grounded ice
%   3 = floating ice
%   4 = Lake Vostok

ocean = bm_mask == 0;
floating = bm_mask == 3;

% Get Thwaites shelf outline
[wx,wy] = antbounds_data('Thwaites','xy');

inside_region = inpolygon(xgrid,ygrid,wx,wy);

% Define the Thwaites floating shelf
thwaites_shelf = ...
    inside_region & ...
    floating & ...
    isfinite(H);

fprintf('Approximate shelf area: %.0f km^2\n', ...
    nnz(thwaites_shelf));

% Create full shelf thickness map
H_shelf_map = H;
H_shelf_map(~thwaites_shelf) = NaN;

%%% Detect the calving front
% A calving-front cell is defined here as:
%   a floating Thwaites shelf cell with at least one neighboring ocean cell.

ocean_neighbors = conv2( ...
    double(ocean), ...
    ones(3), ...
    'same') > 0;

calving_front = thwaites_shelf & ocean_neighbors;

fprintf('Number of detected calving-front cells: %d\n', ...
    nnz(calving_front));

if ~any(calving_front(:))
    error(['No calving-front cells were detected. ', ...
        'Check the Thwaites mask and BedMachine mask.']);
end

%%% Distance inland from the calving front
% bwdist returns Euclidean distance in pixel units.
% Multiply by the physical grid spacing to obtain meters.

distance_from_front = ...
    bwdist(calving_front)*grid_spacing_m;

distance_from_front(~thwaites_shelf) = NaN;

% Near-front region
front_band_width_m = front_band_width_km*1000;

near_front = ...
    thwaites_shelf & ...
    distance_from_front <= front_band_width_m;

H_near_front = H(near_front);

%%% outputs
H_shelf = H(thwaites_shelf);
h_shelf_mean = mean(H_shelf,'omitnan');
h_shelf_eff = mean(H_shelf.^3,'omitnan')^(1/3);

h_front_mean = mean(H_near_front,'omitnan');
h_front_eff = mean(H_near_front.^3,'omitnan')^(1/3);
