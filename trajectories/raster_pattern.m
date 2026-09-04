%% Use region of interest from beam bending model to draw raster pattern

%% Load Coordinate points
load('settling_results.mat');

grid_spacing_km = 1;
grid_spacing_m = grid_spacing_km*1000;

% Extract coordinates
[xgrid, ygrid] = psgrid( ...
    'Thwaites Glacier', 700, grid_spacing_km, 'xy');

% Extract BedMachine thickness
H = bedmachine_interp( ...
    'thickness', ...
    xgrid, ...
    ygrid);


% Find coordinates containing valid shelf data
valid_shelf = isfinite(H_shelf_map);

x_valid = xgrid(valid_shelf);
y_valid = ygrid(valid_shelf);

%% Determine swath range
beta_el = 40; % elevation beam width
beta_az = 40; % azimuth beam width
theta = 40; % look angle
z = 18e3; % altitude of flight

S_range = (z*deg2rad(beta_el)/(cosd(theta)^2)); % swath range

%% Creating Polygons for ROI using Maxwell

% Maximum transect length
max_distance = 20;   % km

% Force everything to column vectors
xg_all = x_gl_match(:);
yg_all = y_gl_match(:);

x_settle_all = maxwell_x_settle(:);
y_settle_all = maxwell_y_settle(:);

% Distance from GL point to Maxwell point
transect_distance = hypot( ...
    x_settle_all - xg_all, ...
    y_settle_all - yg_all) / 1e3;      % km

% Keep only valid transects
valid = isfinite(xg_all) & ...
        isfinite(yg_all) & ...
        isfinite(x_settle_all) & ...
        isfinite(y_settle_all);

xg = xg_all(valid);
yg = yg_all(valid);

xs = x_settle_all(valid);
ys = y_settle_all(valid);

% Normalize strike direction
strike = strike(:) / norm(strike);

% Position of each transect along strike
s = xg * strike(1) + yg * strike(2);

% Sort transects along strike
[~, idx] = sort(s);

xg = xg(idx);
yg = yg(idx);

xs = xs(idx);
ys = ys(idx);

%% Split ROI into two distinct polygons

% Distance between consecutive grounding-line points
d_gl = hypot(diff(xg), diff(yg));

% Find largest gap
[~, split_idx] = max(d_gl);

% Polygon 1
xg_1 = xg(1:split_idx);
yg_1 = yg(1:split_idx);

xs_1 = xs(1:split_idx);
ys_1 = ys(1:split_idx);

x_roi_1 = [xg_1; flipud(xs_1)];
y_roi_1 = [yg_1; flipud(ys_1)];

roi_1 = polyshape(x_roi_1, y_roi_1);


% Polygon 2
xg_2 = xg(split_idx+1:end);
yg_2 = yg(split_idx+1:end);

xs_2 = xs(split_idx+1:end);
ys_2 = ys(split_idx+1:end);

x_roi_2 = [xg_2; flipud(xs_2)];
y_roi_2 = [yg_2; flipud(ys_2)];

roi_2 = polyshape(x_roi_2, y_roi_2);

% Plot
figure;

% Ice shelf thickness
pcolor(xgrid, ygrid, H_shelf_map);
shading flat;

cb = colorbar;
cb.Label.String = 'Ice thickness (m)';

hold on;

% Polygon boundary
[xb1, yb1] = boundary(roi_1);

plot(xb1, yb1, ...
    'k-', ...
    'LineWidth',2, ...
    'DisplayName','Survey ROI');

[xb2, yb2] = boundary(roi_2);

plot(xb2, yb2, ...
    'k-', ...
    'LineWidth',2, ...
    'DisplayName','Survey ROI');


% Grounding line points
scatter(xg, yg, ...
    8, 'k', 'filled', ...
    'DisplayName','Grounding Line');

% Maxwell settling points
scatter(xs, ys, ...
    12, 'r', 'filled', ...
    'DisplayName','Maxwell Settling');
axis equal;
% Add small padding
padding = 5e3;   % km

xlim([min(x_valid) - padding, ...
      max(x_valid) + padding]);

ylim([min(y_valid) - padding, ...
      max(y_valid) + padding]);

%% Create grid tiles in the PCA strike direction and normal direction

% tile resolution
tile_resolution = S_range; 

% Make sure strike is normalized
strike = strike(:) / norm(strike);

% Normal direction: 90 degree rotation
normal = [-strike(2); strike(1)];

% Collect ROI boundary coordinates (ROI 2)

x_all = [xb1(:); xb2(:)];
y_all = [yb1(:); yb2(:)];

% Remove NaNs from polyshape boundaries
valid_boundary = isfinite(x_all) & isfinite(y_all);

x_all = x_all(valid_boundary);
y_all = y_all(valid_boundary);

% Convert ROI coordinates into strike-normal coordinates

% Choose an origin to keep numbers manageable
x0 = mean(x_all);
y0 = mean(y_all);

dx = x_all - x0;
dy = y_all - y0;

% Along-strike coordinate
s_coord = dx * strike(1) + dy * strike(2);

% Cross-strike coordinate
n_coord = dx * normal(1) + dy * normal(2);

%% Create tile edges

s_edges = floor(min(s_coord)/tile_resolution) * tile_resolution : ...
          tile_resolution : ...
          ceil(max(s_coord)/tile_resolution) * tile_resolution;

n_edges = floor(min(n_coord)/tile_resolution) * tile_resolution : ...
          tile_resolution : ...
          ceil(max(n_coord)/tile_resolution) * tile_resolution;

%% Generate tiles

tiles = polyshape.empty;
tile_counter = 0;
roi_total = union(roi_1, roi_2);

occupied = false(length(s_edges)-1, length(n_edges)-1);

for i = 1:length(s_edges)-1

    for j = 1:length(n_edges)-1

        % Tile corners in strike-normal coordinates
        s_tile = [ ...
            s_edges(i), ...
            s_edges(i+1), ...
            s_edges(i+1), ...
            s_edges(i)];

        n_tile = [ ...
            n_edges(j), ...
            n_edges(j), ...
            n_edges(j+1), ...
            n_edges(j+1)];

        % Convert back into x-y coordinates
        x_tile = x0 + ...
            s_tile * strike(1) + ...
            n_tile * normal(1);

        y_tile = y0 + ...
            s_tile * strike(2) + ...
            n_tile * normal(2);

        tile_poly = polyshape(x_tile, y_tile);

        % Check whether tile intersects either ROI
        tile_intersection = intersect(tile_poly, roi_total);
        occupied(i,j) = area(tile_intersection) > 0;

        tile_roi_1 = intersect(tile_poly, roi_1);
        tile_roi_2 = intersect(tile_poly, roi_2);

        if area(tile_roi_1) || area(tile_roi_2) > 0

            tile_counter = tile_counter + 1;
            tiles(tile_counter) = tile_poly;


        end
    end
end

fprintf('Number of tiles = %d\n', tile_counter);

%% Plot tiles
figure;
hold on;

% Ice shelf thickness
pcolor(xgrid, ygrid, H_shelf_map);
shading flat;

cb = colorbar;
cb.Label.String = 'Ice thickness (m)';

% Polygon boundary
[xb1, yb1] = boundary(roi_1);

plot(xb1, yb1, ...
    'k-', ...
    'LineWidth',2, ...
    'DisplayName','Survey ROI');

[xb2, yb2] = boundary(roi_2);

plot(xb2, yb2, ...
    'k-', ...
    'LineWidth',2, ...
    'DisplayName','Survey ROI');


% Grounding line points
scatter(xg, yg, ...
    8, 'k', 'filled', ...
    'DisplayName','Grounding Line');

% Maxwell settling points
scatter(xs, ys, ...
    12, 'r', 'filled', ...
    'DisplayName','Maxwell Settling');


for k = 1:length(tiles)

    plot(tiles(k), ...
        'FaceColor','none', ...
        'EdgeColor','r', ...
        'LineWidth',2);

end
axis equal;

% Add small padding
padding = 30e3;   % km

xlim([min(x_valid) - padding, ...
      max(x_valid) + padding]);

ylim([min(y_valid) - padding, ...
      max(y_valid) + padding]);


%% ============================================================
%  RASTER PATH + TILE SURVEY TIMES
% ============================================================

V_ground = 20;       % aircraft ground speed [m/s]
dt = 5;              % trajectory time step [s]

raster_mode = 'strike';   % 'strike' or 'normal'

% One survey time for every occupied tile
tile_times = nan(size(occupied));    % [hr]

% Track total distance traveled while constructing path
distance_traveled = 0;               % [m]


%% ============================================================
%  STRIKE-DIRECTION RASTER
% ============================================================

if strcmpi(raster_mode,'strike')

    s_path = [];
    n_path = [];

    direction = 1;

    for j = 1:size(occupied,2)

        % Tiles that exist in this normal-direction row
        active_i = find(occupied(:,j));

        if isempty(active_i)
            continue
        end

        % Ends of raster leg
        s_min = s_edges(min(active_i));
        s_max = s_edges(max(active_i)+1);

        % Aircraft flies on upper normal edge of row
        n_current = n_edges(j+1);

        % Determine direction of this pass
        if direction == 1
            s_start = s_min;
            s_end   = s_max;
        else
            s_start = s_max;
            s_end   = s_min;
        end


        % Move aircraft to beginning of survey leg

        if isempty(s_path)

            % First point of entire flight
            s_path = s_start;
            n_path = n_current;

        else

            s_prev = s_path(end);
            n_prev = n_path(end);

            % First move along strike until aligned
            if abs(s_prev - s_start) > 1e-6

                distance_traveled = distance_traveled + ...
                    abs(s_start - s_prev);

                s_path(end+1) = s_start;
                n_path(end+1) = n_prev;

            end

            % Then move perpendicular to strike
            if abs(n_prev - n_current) > 1e-6

                distance_traveled = distance_traveled + ...
                    abs(n_current - n_prev);

                s_path(end+1) = s_start;
                n_path(end+1) = n_current;

            end

        end


        % Assign survey time to each tile in this row

        for ii = active_i(:)'

            % Center of this tile along strike
            s_center = ...
                (s_edges(ii) + s_edges(ii+1))/2;

            % Distance from beginning of this pass
            distance_into_leg = ...
                abs(s_center - s_start);

            % Total distance traveled when aircraft reaches tile center
            distance_at_tile = ...
                distance_traveled + distance_into_leg;

            % Survey time [hr]
            tile_times(ii,j) = ...
                distance_at_tile / V_ground / 3600;

        end


        % Fly survey leg

        survey_distance = abs(s_end - s_start);

        distance_traveled = ...
            distance_traveled + survey_distance;

        s_path(end+1) = s_end;
        n_path(end+1) = n_current;


        % Reverse direction for next row
        direction = -direction;

    end


% ============================================================
%  NORMAL-DIRECTION RASTER
% ============================================================

elseif strcmpi(raster_mode,'normal')

    s_path = [];
    n_path = [];

    direction = 1;

    for i = 1:size(occupied,1)

        % Tiles that exist in this strike-direction column
        active_j = find(occupied(i,:));

        if isempty(active_j)
            continue
        end

        % Ends of raster leg
        n_min = n_edges(min(active_j));
        n_max = n_edges(max(active_j)+1);

        % Aircraft flies along one strike-direction tile edge
        s_current = s_edges(i);

        % Determine pass direction
        if direction == 1
            n_start = n_min;
            n_end   = n_max;
        else
            n_start = n_max;
            n_end   = n_min;
        end


        % Move aircraft to beginning of survey leg

        if isempty(s_path)

            % First point of entire flight
            s_path = s_current;
            n_path = n_start;

        else

            s_prev = s_path(end);
            n_prev = n_path(end);

            % First move along normal until aligned
            if abs(n_prev - n_start) > 1e-6

                distance_traveled = distance_traveled + ...
                    abs(n_start - n_prev);

                s_path(end+1) = s_prev;
                n_path(end+1) = n_start;

            end

            % Then move along strike to next column
            if abs(s_prev - s_current) > 1e-6

                distance_traveled = distance_traveled + ...
                    abs(s_current - s_prev);

                s_path(end+1) = s_current;
                n_path(end+1) = n_start;

            end

        end


        % Assign survey time to each tile in this column

        for jj = active_j(:)'

            % Center of tile along normal direction
            n_center = ...
                (n_edges(jj) + n_edges(jj+1))/2;

            % Distance from beginning of this pass
            distance_into_leg = ...
                abs(n_center - n_start);

            % Total distance traveled
            distance_at_tile = ...
                distance_traveled + distance_into_leg;

            % Survey time [hr]
            tile_times(i,jj) = ...
                distance_at_tile / V_ground / 3600;

        end


        % Fly survey leg

        survey_distance = abs(n_end - n_start);

        distance_traveled = ...
            distance_traveled + survey_distance;

        s_path(end+1) = s_current;
        n_path(end+1) = n_end;


        % Reverse direction
        direction = -direction;

    end

else

    error('raster_mode must be ''strike'' or ''normal''.');

end


%% ============================================================
%  CONVERT s-n PATH TO x-y
% ============================================================

x_path = x0 + ...
    s_path*strike(1) + ...
    n_path*normal(1);

y_path = y0 + ...
    s_path*strike(2) + ...
    n_path*normal(2);

x_path = x_path(:);
y_path = y_path(:);


%% ============================================================
%  PARAMETERIZE FLIGHT PATH IN TIME
% ============================================================

dx = diff(x_path);
dy = diff(y_path);

segment_length = hypot(dx,dy);

distance_path = [0; cumsum(segment_length)];

total_distance = distance_path(end);

fprintf('Total flight distance = %.2f km\n', ...
    total_distance/1e3);


% Time at raster vertices
time_path = distance_path / V_ground;

total_time = time_path(end);

fprintf('Total flight time = %.2f hr\n', ...
    total_time/3600);


% Fine time vector
time_flight = (0:dt:total_time)';

if time_flight(end) < total_time
    time_flight(end+1) = total_time;
end


% Interpolated aircraft coordinates
x_flight = interp1( ...
    time_path, ...
    x_path, ...
    time_flight, ...
    'linear');

y_flight = interp1( ...
    time_path, ...
    y_path, ...
    time_flight, ...
    'linear');


%% ============================================================
%  PLOT
% ============================================================

figure;
hold on;

% Shelf background

shelf_plot = ones(size(H_shelf_map));
shelf_plot(~isfinite(H_shelf_map)) = NaN;

h = pcolor(xgrid, ygrid, shelf_plot);
shading flat;

gray_color = [0.8 0.8 0.8];

rgb = repmat( ...
    reshape(gray_color,1,1,3), ...
    size(shelf_plot,1), ...
    size(shelf_plot,2), ...
    1);

set(h,'CData',rgb);

set(h, ...
    'AlphaData',double(isfinite(H_shelf_map)), ...
    'FaceAlpha','flat');


% ROI outlines

plot(roi_1, ...
    'FaceColor','none', ...
    'EdgeColor','k', ...
    'LineWidth',2);

plot(roi_2, ...
    'FaceColor','none', ...
    'EdgeColor','k', ...
    'LineWidth',2);


% ============================================================
%  TILE COLORS
% ============================================================

cmap = parula(256);
colormap(cmap);

% Only occupied/surveyed tiles
valid_times = tile_times(isfinite(tile_times));

tmin = min(valid_times);
tmax = max(valid_times);

% Make colorbar use full available range
clim([tmin tmax]);


for i = 1:size(occupied,1)

    for j = 1:size(occupied,2)

        if ~occupied(i,j)
            continue
        end

        if isnan(tile_times(i,j))
            continue
        end


        % Construct this tile directly from its grid edges

        s_tile = [ ...
            s_edges(i), ...
            s_edges(i+1), ...
            s_edges(i+1), ...
            s_edges(i)];

        n_tile = [ ...
            n_edges(j), ...
            n_edges(j), ...
            n_edges(j+1), ...
            n_edges(j+1)];


        % Convert tile to x-y
        x_tile = x0 + ...
            s_tile*strike(1) + ...
            n_tile*normal(1);

        y_tile = y0 + ...
            s_tile*strike(2) + ...
            n_tile*normal(2);


        % Convert tile time into colormap color

        cval = ...
            (tile_times(i,j) - tmin) / ...
            (tmax - tmin);

        cval = max(0,min(1,cval));

        color_idx = ...
            round(cval*(size(cmap,1)-1)) + 1;

        tile_color = cmap(color_idx,:);


        % Draw tile

        tile_poly = polyshape(x_tile,y_tile);

        plot(tile_poly, ...
            'FaceColor',tile_color, ...
            'EdgeColor',[0.6 0.6 0.6], ...
            'LineWidth',0.5, ...
            'HandleVisibility','off');

    end

end


% Flight path

plot(x_path, y_path, ...
    'r-', ...
    'LineWidth',1.5, ...
    'HandleVisibility','off');


% Start point

scatter(x_flight(1), y_flight(1), ...
    80, ...
    'g', ...
    'filled', ...
    'MarkerEdgeColor','k', ...
    'DisplayName','Start');


% End point
scatter(x_flight(end), y_flight(end), ...
    80, ...
    'r', ...
    'filled', ...
    'MarkerEdgeColor','k', ...
    'DisplayName','End');

% Add flight-direction arrow
arrow_fracs = [0.2 0.4 0.6 0.8];

arrow_length = 4e3;   % arrowhead length [m]
arrow_width  = 2e3;   % arrowhead width [m]

for f = arrow_fracs

    arrow_idx = round(length(x_flight) * f);

    % Direction of flight
    dx = x_flight(arrow_idx+1) - x_flight(arrow_idx);
    dy = y_flight(arrow_idx+1) - y_flight(arrow_idx);

    % Normalize direction vector
    L = hypot(dx,dy);

    ux = dx/L;
    uy = dy/L;

    % Perpendicular direction
    px = -uy;
    py = ux;

    % Location of arrow tip
    xtip = x_flight(arrow_idx);
    ytip = y_flight(arrow_idx);

    % Back center of triangle
    xback = xtip - arrow_length*ux;
    yback = ytip - arrow_length*uy;

    % Triangle corners
    x_arrow = [ ...
        xtip, ...
        xback + arrow_width/2*px, ...
        xback - arrow_width/2*px];

    y_arrow = [ ...
        ytip, ...
        yback + arrow_width/2*py, ...
        yback - arrow_width/2*py];

    % Draw arrowhead
    patch(x_arrow, y_arrow, 'k', ...
        'EdgeColor','k', ...
        'HandleVisibility','off');

end

% Colorbar

cb = colorbar;
cb.Label.String = 'Tile survey time (hr)';


% Figure formatting

axis equal;

xlabel('x (m)');
ylabel('y (m)');

title(sprintf( ...
    '%s Raster, Total Flight Time %.2f hr', ...
    upper(raster_mode), ...
    total_time/3600));


padding = 30e3;   % [m]

xlim([ ...
    min(x_valid)-padding, ...
    max(x_valid)+padding]);

ylim([ ...
    min(y_valid)-padding, ...
    max(y_valid)+padding]);


figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/raster_pattern/20km';

% Create folder if it doesn't exist
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

exportgraphics(gcf, ...
    fullfile(figure_dir, 'colored_tiles_strike.png'), ...
    'Resolution', 300);

% exportgraphics(gcf, ...
%    fullfile(figure_dir, 'colored_tiles_normal.png'), ...
%    'Resolution', 300);
%% Battery simulation along flight path
battery_capacity = 27.2;      % kWh
E_initial = battery_capacity; % start fully charged

dt = time_flight(2) - time_flight(1);   % seconds

% Preallocate
P_solar_flight = zeros(size(time_flight));
P_usage_flight = zeros(size(time_flight));
P_net_flight   = zeros(size(time_flight));

E_battery = zeros(size(time_flight));
E_battery(1) = E_initial;

% Absolute flight time

% Example mission start time
date_start = datenum(2026,1,20,0,0,0);

% Convert elapsed seconds -> MATLAB datenum
time_num_flight = date_start + time_flight/(24*3600);

% March through flight
for k = 1:length(time_flight)-1

    % -------------------------------------------------
    % Solar power at current aircraft position/time
    % -------------------------------------------------

    % If x_flight/y_flight are Antarctic polar stereographic,
    % convert back to lat/lon
    [lat_k, lon_k] = ps2ll(x_flight(k), y_flight(k));

    P_solar_flight(k) = calculate_solar_power( ...
        lat_k, ...
        lon_k, ...
        time_num_flight(k));

    % -------------------------------------------------
    % Aircraft power usage
    % -------------------------------------------------

    V_ground = 20;  % m/s
    V_wind   = 0;   % m/s

    P_usage_flight(k) = calculate_flight_power( ...
        V_ground, ...
        V_wind);

    % -------------------------------------------------
    % Net battery power
    % -------------------------------------------------

    P_net_flight(k) = ...
        P_solar_flight(k) - P_usage_flight(k);

    % -------------------------------------------------
    % Update battery energy
    % -------------------------------------------------

    dt_k = time_flight(k+1) - time_flight(k);

    E_battery(k+1) = E_battery(k) + ...
        P_net_flight(k)*dt_k/(1000*3600);

    % Battery cannot exceed capacity
    E_battery(k+1) = min( ...
        E_battery(k+1), ...
        battery_capacity);

    % Battery cannot go below zero
    E_battery(k+1) = max( ...
        E_battery(k+1), ...
        0);
end

% convert into percentage
E_percent = 100* E_battery / battery_capacity; 
%% generate plots

figure;
hold on;

% plot the shelf as just an outline
shelf_plot = ones(size(H_shelf_map));
shelf_plot(~isfinite(H_shelf_map)) = NaN;

h = pcolor(xgrid, ygrid, shelf_plot);
shading flat;

gray_color = [0.8 0.8 0.8];
% Convert CData to truecolor RGB so it does not use the colormap
rgb = repmat(reshape(gray_color,1,1,3), ...
             size(shelf_plot,1), ...
             size(shelf_plot,2), ...
             1);
set(h, 'CData', rgb);
% Make NaN regions transparent
set(h, ...
    'AlphaData', double(isfinite(H_shelf_map)), ...
    'FaceAlpha','flat');

% ROI
plot(roi_1, ...
    'FaceColor','none', ...
    'EdgeColor','k', ...
    'LineWidth',2);

plot(roi_2, ...
    'FaceColor','none', ...
    'EdgeColor','k', ...
    'LineWidth',2);

% Tiles
for k = 1:length(tiles)

    plot(tiles(k), ...
        'FaceColor','none', ...
        'EdgeColor',[0.7 0.7 0.7], ...
        'LineWidth',2);

end


plot(x_path, y_path, ...
    'k-', ...
    'LineWidth',2);

scatter(x_flight, y_flight, ...
    5, ...
    E_percent, ...
    'filled');

% Start point
scatter(x_flight(1), y_flight(1), ...
    80, ...
    'g', ...
    'filled', ...
    'MarkerEdgeColor','k', ...
    'DisplayName','Start');

% End point
scatter(x_flight(end), y_flight(end), ...
    80, ...
    'r', ...
    'filled', ...
    'MarkerEdgeColor','k', ...
    'DisplayName','End');

% Add flight-direction arrow
arrow_fracs = [0.2 0.4 0.6 0.8];

arrow_length = 4e3;   % arrowhead length [m]
arrow_width  = 2e3;   % arrowhead width [m]

for f = arrow_fracs

    arrow_idx = round(length(x_flight) * f);

    % Direction of flight
    dx = x_flight(arrow_idx+1) - x_flight(arrow_idx);
    dy = y_flight(arrow_idx+1) - y_flight(arrow_idx);

    % Normalize direction vector
    L = hypot(dx,dy);

    ux = dx/L;
    uy = dy/L;

    % Perpendicular direction
    px = -uy;
    py = ux;

    % Location of arrow tip
    xtip = x_flight(arrow_idx);
    ytip = y_flight(arrow_idx);

    % Back center of triangle
    xback = xtip - arrow_length*ux;
    yback = ytip - arrow_length*uy;

    % Triangle corners
    x_arrow = [ ...
        xtip, ...
        xback + arrow_width/2*px, ...
        xback - arrow_width/2*px];

    y_arrow = [ ...
        ytip, ...
        yback + arrow_width/2*py, ...
        yback - arrow_width/2*py];

    % Draw arrowhead
    patch(x_arrow, y_arrow, 'k', ...
        'EdgeColor','k', ...
        'HandleVisibility','off');

end

axis equal;

colormap(flipud(nebula));
cb = colorbar;
cb.Label.String = 'Battery percent';

xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Battery Life Raster Pattern, Total flight time %.2f hours', total_time/3600));

% Add small padding
padding = 30e3;   % km

xlim([min(x_valid) - padding, ...
      max(x_valid) + padding]);

ylim([min(y_valid) - padding, ...
      max(y_valid) + padding]);

% exportgraphics(gcf, ...
%     fullfile(figure_dir, 'strike_battery.png'), ...
%     'Resolution', 300);

% exportgraphics(gcf, ...
%    fullfile(figure_dir, 'normal_battery.png'), ...
%    'Resolution', 300);

%% Function for calculating solar flux 
function [P_solar, Id, Is, zeta, solar_az] = calculate_solar_power(lat, lon, time_num)
% CALCULATE_SOLAR_POWER Calculate solar power available to HAPS
%
% Inputs:
%   lat      - latitude [deg]
%   lon      - longitude [deg]
%   time_num - MATLAB datenum in UTC
%
% Outputs:
%   P_solar  - electrical solar power [W]
%   Id       - direct irradiance on panel [W/m^2]
%   Is       - diffuse irradiance on panel [W/m^2]
%   zeta     - solar zenith angle [deg]
%   solar_az - solar azimuth angle [deg]


% Solar position

time_zone = 0;   % UTC
rotation  = 0;   % coordinate-system rotation [deg]
dst       = false;

[angles, ~] = solarPosition( ...
    time_num, ...
    lat, ...
    lon, ...
    time_zone, ...
    rotation, ...
    dst);

zeta     = angles(:,1);   % solar zenith angle [deg]
solar_az = angles(:,2);   % solar azimuth [deg]


% Extraterrestrial irradiance

% Convert datenum to datetime
date = datetime(time_num, 'ConvertFrom', 'datenum');

% Day of year
n = day(date, 'dayofyear');

% Solar constant (Spencer 1971)
I_sc = 1361;   % W/m^2

% Earth-Sun distance correction
Gamma = 2*pi*(n - 1)/365;

E0 = ...
    1.000110 + ...
    0.034221*cos(Gamma) + ...
    0.001280*sin(Gamma) + ...
    0.000719*cos(2*Gamma) + ...
    0.000077*sin(2*Gamma);

% Extraterrestrial normal irradiance
I0_n = I_sc .* E0;   % W/m^2

% Solar panel parameters
panel_area = 33.2;    % solar panel area [m^2] (0.8*wing_area)
eta_panel  = 0.216;  % panel efficiency [-] (Dewald 2024)
eta_MTTP = 0.9615;
eta_wiring = 0.99;
eta_gearbox = 0.986;

eta_total = eta_panel*eta_MTTP*eta_wiring*eta_gearbox;

% Atmospheric transmittances (find source)
Tr = 0.99;   % Rayleigh scattering transmittance [-]
To = 0.97;   % ozone transmittance [-]
Tg = 0.99;   % mixed-gas transmittance [-]
T_atm = Tr*To*Tg; 

% Panel orientation
beta  = 0;   % panel tilt from horizontal [deg]
gamma = 0;   % panel azimuth [deg]

alpha = solar_az;

% Cosine of solar incidence angle
cos_theta = ...
    cosd(zeta).*cosd(beta) + ...
    sind(zeta).*sind(beta).*cosd(alpha - gamma);

% No direct sunlight on backside of panel
cos_theta = max(cos_theta, 0);

% Direct irradiance
Id = I0_n .* T_atm .* cos_theta;

% Diffuse irradiance

% Prevent negative diffuse irradiance when Sun is below horizon
cos_zeta = max(cosd(zeta), 0);
Is = ...
    0.5 .* I0_n .* cos_zeta .* ...
    To .* Tg .* (1 - Tr);

% Available solar power
P_solar = panel_area .* eta_total .* (Id + Is);

end

%% Flight Power consumption
function [P_el, V_air, Cl, Cd, Ft] = calculate_flight_power(V_ground, V_wind)
% CALCULATE_FLIGHT_POWER Calculate electrical power required for flight
%
% Inputs:
%   V_ground - Ground speed [m/s]
%   V_wind   - Wind speed [m/s]
%              positive = tailwind
%              negative = headwind
%
% Outputs:
%   P_el  - Electrical power required [W]
%   V_air - Airspeed [m/s]
%   Cl    - Lift coefficient
%   Cd    - Drag coefficient
%   Ft    - Required thrust [N]

% Aircraft parameters
rho_air   = 0.1216;   % air density [kg/m^3]
wing_area = 41.5;     % wing area [m^2]
eta_prop  = 0.9;      % propeller efficiency
g         = 9.81;     % gravity [m/s^2]

TOGW = 214;            % takeoff gross mass [kg]
LD   = 28.4;           % lift-to-drag ratio

% Airspeed
V_air = V_ground - V_wind;

if V_air <= 0
    error('Airspeed must be greater than zero.');
end

% Lift coefficient
Cl = (2*TOGW*g) / ...
     (rho_air*V_air^2*wing_area);

% Drag coefficient
Cd = Cl / LD;

% Required thrust
Ft = 0.5*rho_air*V_air^2*wing_area*Cd;

% Electrical power
P_el = V_air*Ft/eta_prop;

end