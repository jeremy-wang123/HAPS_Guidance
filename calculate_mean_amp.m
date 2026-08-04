function [mean_amp] = calculate_mean_amp(glacier_name, mode)
%EXTRACT_MEAN_AMP extracts the mean tidal amplitude of a region
%   Inputs: 
%   - Glacier name (i.e. Thwaites)
%   - mode ('top', 'top2', 'all')
%   Output: mean tidal amplitude of the ice shelf region
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008')
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008';

boundary_name = glacier_name;                 % 'Thwaites'
region_name   = [glacier_name ' Glacier'];    % 'Thwaites Glacier'

% extract grids
[xgrid,ygrid] = psgrid(region_name,700,1,'xy');

% extract the amplitudes
[constit_list, constit_mean_amp, constit_idx] = select_constituent_amplitudes(region_name, ...
    Model, mode);

[lat,lon] = ps2ll(xgrid,ygrid);
[amps,~,~,~] = tmd_extract_HC(Model,lat,lon,'z',constit_idx);

amp_sum = squeeze(sum(amps,1));

% Get Thwaites boundary
[wx,wy] = antbounds_data(boundary_name,'xy');

% Points inside the Thwaites polygon
inside_region = inpolygon(xgrid,ygrid,wx,wy);

% Identify floating ice
shelf = isiceshelf(xgrid,ygrid);

% Thwaites floating shelf mask
thwaites_shelf = inside_region & shelf;

x_shelf = xgrid;
y_shelf = ygrid;

x_shelf(~thwaites_shelf) = NaN;
y_shelf(~thwaites_shelf) = NaN;

% Calculate the mean amplitude for the floating ice shelf
mean_amp = mean(amp_sum(thwaites_shelf), 'omitnan');

