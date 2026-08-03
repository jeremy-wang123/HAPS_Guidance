%% Analytical solution to elastic beam bending

%% Parameter values 
%{
wa = tidal amplitude at the edge of the ice shelf (determined from the
CATS2008 model
mu = Poisson's ratio (homogenous for ice)
E = Young's Modulus (uniform for ice)
I = moment of inertia of the section per unit width (see note)
rho = density of water
g = acceleration from gravity

For our purposes, we can treat I  = ∫ y^2 dy from -h/2 to h/2 = h^3 / 12
where h is the thickness of the ice beam. This is assumes uniform thickness
of a rectangular ice shelf and symmetry about y=0

This approximation works for only long shelves (where L >>
5π/4(lambda)

For ice shelves less than that length, we need to treat them as
elastic-plastic deformation
%}

%% determine ice shelf thickness using mean ice shelf thickness
% Morlighem et al., 2017)

addpath('/Users/jeremywang/Documents/MATLAB/BedMachine');
[xgrid,ygrid] = psgrid('Thwaites Glacier',700,1,'xy');
H = bedmachine_interp('thickness',xgrid,ygrid);

% Get Thwaites boundary
[wx,wy] = antbounds_data('Thwaites','xy');

% Points inside the Thwaites polygon
inside_region = inpolygon(xgrid,ygrid,wx,wy);

% Identify floating ice
shelf = isiceshelf(xgrid,ygrid);

% Thwaites floating shelf mask
thwaites_shelf = inside_region & shelf & isfinite(H);

% plotting
H_shelf_map = H;
H_shelf_map(~thwaites_shelf) = NaN;

figure
pcolor(xgrid,ygrid,H_shelf_map)
shading flat
axis image
hold on

% plot(wx,wy,'k','LineWidth',2)     % Thwaites boundary
antbounds('gl','k','LineWidth',1) % Grounding line (optional)

colormap(parula)
cb = colorbar;
ylabel(cb,'Thickness (m)')
xlabel('x (m)')
ylabel('y (m)')
title('BedMachine Thickness of the Thwaites Ice Shelf')

H_shelf = H(thwaites_shelf); % Extract thickness values for the floating shelf
% for thickness profiles, mean, medium, and effective thickness calculated
h_mean   = mean(H_shelf,'omitnan');
h_median = median(H_shelf,'omitnan');
h_eff    = mean(H_shelf.^3,'omitnan')^(1/3); % take the mean of h^3 instead of h

%% Thwaites Ice Shelf thickness profile from the calving front
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

%%% summary statistics
H_shelf = H(thwaites_shelf);

h_shelf_mean = mean(H_shelf,'omitnan');
h_shelf_median = median(H_shelf,'omitnan');
h_shelf_eff = mean(H_shelf.^3,'omitnan')^(1/3);

h_front_mean = mean(H_near_front,'omitnan');
h_front_median = median(H_near_front,'omitnan');
h_front_eff = mean(H_near_front.^3,'omitnan')^(1/3);

fprintf('\nWhole Thwaites shelf:\n');
fprintf('Mean thickness:                 %.1f m\n',h_shelf_mean);
fprintf('Median thickness:               %.1f m\n',h_shelf_median);
fprintf('Rigidity-equivalent thickness:  %.1f m\n',h_shelf_eff);

fprintf('\nWithin %.0f km of calving front:\n', ...
    front_band_width_km);

fprintf('Mean thickness:                 %.1f m\n',h_front_mean);
fprintf('Median thickness:               %.1f m\n',h_front_median);
fprintf('Rigidity-equivalent thickness:  %.1f m\n',h_front_eff);

%% Determine wa with CATS model
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008')
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008';

% extract grids
[xgrid,ygrid] = psgrid('Thwaites Glacier',700,1,'xy');

% extract the amplitudes
[constit_list, constit_mean_amp, constit_idx] = select_constituent_amplitudes('Thwaites Glacier', ...
    Model, 'top2');

[lat,lon] = ps2ll(xgrid,ygrid);
[amps,~,~,~] = tmd_extract_HC(Model,lat,lon,'z',constit_idx);

amp_sum = squeeze(sum(amps,1));

% Get Thwaites boundary
[wx,wy] = antbounds_data('Thwaites','xy');

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

%% Solving w(x)

% Define parameters
mu = 0.3; % poisson's ratio
rho = 1028; % kg/m^3
g = 9.81; % m/s^2
E = 2.7e4 * 10^5; % Pa
h = h_front_eff; % use the mean effective thickness of calving front
I = h^3/12;
wa = mean_amp;

lambda = (rho*g*(1-mu^2)/(4*E*I))^(1/4);
min_L = 5*pi/(4*lambda); % minimum ice shelf length for elastic estimation to make sense (according to Holdsworth)

spacing = 100;
L = 15e3; % length of ice shelf 
x = 0:spacing:L; 
x_km = x/1000;

w = wa * (1 - (exp(-lambda*x) .* (cos(lambda*x) + sin(lambda*x))));

%% Settling-distance criterion based on displacement uncertainty

% Vertical displacement uncertainty
% Must have the same units as w and wa, normally meters
sigma_z = 0.002;  % m

% calculate settling distance analytically
x_settle = log(sqrt(2)*abs(wa)/sigma_z)/lambda;
x_settle_km = x_settle/1000;

fprintf('1/lambda flexural length: %.2f km\n', ...
    1/lambda/1000);

fprintf('Displacement uncertainty sigma_z: %.4f m\n', ...
    sigma_z);

fprintf('Settling distance within +/- sigma_z: %.2f km\n', ...
    x_settle_km);

%% Plotting profile

figure;

plot(x_km, w, 'LineWidth', 2);
hold on;

% Far-field steady-state tidal displacement
yline(wa, '--', 'Far-field tide', ...
    'LabelHorizontalAlignment', 'left');

% Upper and lower uncertainty bounds around steady state
upper_bound = wa + sigma_z;
lower_bound = wa - sigma_z;

yline(upper_bound, ':', ...
    '$\pm\sigma_z$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right');

yline(lower_bound, ':');

% Settling-distance marker
xline(x_settle_km, ':', ...
    sprintf('Settling distance = %.1f km', x_settle_km), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Elastic tidal flexure of the Thwaites Ice Shelf');

xlim([0, L/1000]);

% Include the displacement profile and uncertainty bounds in the y-limits
all_y_values = [w(:); upper_bound; lower_bound];

w_range = max(all_y_values) - min(all_y_values);

% Avoid zero padding if the plotted values are nearly constant
if w_range == 0
    vertical_padding = max(abs(all_y_values))*0.1;
else
    vertical_padding = 0.10*w_range;
end

ylim([min(all_y_values) - vertical_padding, ...
      max(all_y_values) + vertical_padding]);

grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);

% figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures';
% exportgraphics(gcf, fullfile(figure_dir,'elastic_beam_profile.jpg'), ...
% 'Resolution',300);

%%
function [selected_names, selected_amps, selected_idx] = ...
    select_constituent_amplitudes(glacier_name, Model, mode, specific_constituent)

% SELECT_CONSTITUENT_AMPLITUDES Select tidal constituents by amplitude.
%
% Inputs:
%   glacier_name          Glacier name used by extract_mean_amp
%   Model                 Tidal model path or model identifier
%   mode                  'top', 'top2', 'all', or 'specific'
%   specific_constituent  Constituent name for 'specific' mode
%
% Outputs:
%   selected_names        Cell array of selected constituent names
%   selected_amps         Corresponding mean amplitudes
%   selected_idx          Indices in conlist
%
% Examples:
%   select_constituent_amplitudes( ...
%       'Thwaites Glacier', Model, 'top')
%
%   select_constituent_amplitudes( ...
%       'Thwaites Glacier', Model, 'top2')
%
%   select_constituent_amplitudes( ...
%       'Thwaites Glacier', Model, 'all')
%
%   select_constituent_amplitudes( ...
%       'Thwaites Glacier', Model, 'specific', 'K1')

    if nargin < 4
        specific_constituent = '';
    end

    % Get mean amplitudes
    mean_amps = extract_mean_amp(glacier_name);

    names = fieldnames(mean_amps);
    values = cell2mat(struct2cell(mean_amps));

    % Ensure column vectors
    names = names(:);
    values = values(:);

    % Sort amplitudes from largest to smallest
    [sorted_values, sort_idx] = sort(values,'descend');
    sorted_names = names(sort_idx);

    % Select constituents
    switch lower(mode)

        case 'top'
            selected_names = sorted_names(1);
            selected_amps = sorted_values(1);

        case 'top2'
            n_select = min(2,numel(sorted_names));

            selected_names = sorted_names(1:n_select);
            selected_amps = sorted_values(1:n_select);

        case 'all'
            selected_names = sorted_names;
            selected_amps = sorted_values;

        case 'specific'
            if isempty(specific_constituent)
                error(['A constituent name must be provided ', ...
                       'when mode is ''specific''.']);
            end

            match_idx = find(strcmpi(names,specific_constituent),1);

            if isempty(match_idx)
                error('Constituent "%s" was not found.', ...
                    specific_constituent);
            end

            selected_names = names(match_idx);
            selected_amps = values(match_idx);

        otherwise
            error(['Invalid mode "%s". Choose ''top'', ''top2'', ', ...
                   '''all'', or ''specific''.'],mode);
    end

    % Match selected names to model constituent list
    conlist = extract_conlist(Model);

    % extract_conlist may return a character array
    if ischar(conlist)
        conlist = cellstr(conlist);
    end

    conlist = strtrim(conlist);

    selected_idx = NaN(numel(selected_names),1);

    for k = 1:numel(selected_names)

        idx = find(strcmpi(conlist,selected_names{k}),1);

        if isempty(idx)
            warning('Constituent "%s" was not found in conlist.', ...
                selected_names{k});
        else
            selected_idx(k) = idx;
        end
    end
end