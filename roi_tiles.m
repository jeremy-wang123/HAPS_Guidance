%% Tile region of interest
% create grid of 10km seaward and inland at each grounding line location.
clear;clc;

%% extract grounding line coordinates

% adding model
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008');
addpath('/Users/jeremywang/Documents/MATLAB/measures_v3.1.2/measures');
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008'; 

% extract the coordinates of Thwaites glacier
figure;
mapzoomps('Thwaites Glacier');
x_limits = xlim;
y_limits = ylim;
close;

% specify grid resolution (much be smaller than SAR radar range)
grid_resolution = 5e3; % 5 km
radius = 10e3; % 10 km range for grounding zone
search_spacing = 500; 

% create the raster grid
x_grid = x_limits(1):grid_resolution:x_limits(2);
y_grid = y_limits(1):grid_resolution:y_limits(2);
[X_grid, Y_grid] = meshgrid(x_grid, y_grid);

% Grounding line in the broader plotting region
[xgl_all, ygl_all] = measures_data( ...
    'gl', x_limits, y_limits, 'xy');
%% Method 1: Extract entire Region
% iterate through each gridpoint, determine if it switches status any 
% point within 10km radius

sampling_zone = false(size(X_grid));

for idx = 1:numel(X_grid)

    x0 = X_grid(idx);
    y0 = Y_grid(idx);

    sampling_zone(idx) = is_near_grounding_line( ...
        x0, y0, radius, search_spacing);
end

figure;

% Plot a constant-color layer
shade_layer = ones(size(sampling_zone));

h = imagesc(x_grid, y_grid, shade_layer);
set(gca, 'YDir', 'normal');

% Make only sampling_zone == true visible
h.AlphaData = 0.35 * double(sampling_zone);
h.AlphaDataMapping = 'none';

% Set overlay color
colormap(gca, [1 0 0]);

% Plot grounding line above overlay
hold on;
measuresps('gl', 'k', 'LineWidth', 1.5);

xlim(x_limits);
ylim(y_limits);
xlabel('x (m)');
ylabel('y (m)');
axis image;

%% Method 2: Extract Individual Glacier

% drainage basin
[xbasin, ybasin] = basin_data( ...
    'imbie refined', 'Thwaites', 'xy');

% Keep grounding-line points inside basin
inside = inpolygon( ...
    xgl_all, ygl_all, xbasin, ybasin);

xgl_glacier = xgl_all;
ygl_glacier = ygl_all;

% Preserve separation between line segments
xgl_glacier(~inside) = NaN;
ygl_glacier(~inside) = NaN;

% Find distance from each grid point to Thwaites grounding line
grid_points = [X_grid(:), Y_grid(:)];
gl_points   = [xgl_glacier(:), ygl_glacier(:)];

% Distance to nearest grounding-line point
[~, distance_to_gl] = knnsearch(gl_points, grid_points);

% Logical mask: true if within radius
sampling_zone = reshape( ...
    distance_to_gl <= radius, ...
    size(X_grid));

figure;
mapzoomps('Thwaites Glacier')
measuresps('gl', 'k')

hold on;

% Plot a constant-color layer
shade_layer = ones(size(sampling_zone));

h = imagesc(x_grid, y_grid, shade_layer);
set(gca, 'YDir', 'normal');

% Make only sampling_zone == true visible
h.AlphaData = 0.35 * double(sampling_zone);
h.AlphaDataMapping = 'none';

% Set overlay color
colormap(gca, [1 0 0]);

% Plot grounding line above overlay
plot(xgl_glacier, ygl_glacier, 'k-', 'LineWidth', 3)

xlim(x_limits);
ylim(y_limits);
xlabel('x (m)');
ylabel('y (m)');
axis image;

%% PCA for singular flight path

% establish the swath range
beta_el = 40; % elevation beam width
beta_az = 40; % azimuth beam width
theta = 40; % look angle
z = 18e3; % altitude of flight

S_range = (z*deg2rad(beta_el)/(cosd(theta)^2)); % in m

% Extract actual x and y coordinates in meters
x_gl = X_grid(sampling_zone);
y_gl = Y_grid(sampling_zone);

% conduct pca
P = [x_gl(:), y_gl(:)];
P = P(all(isfinite(P), 2), :);
[coeff, score, latent, ~, explained, mu] = pca(P);

direction = coeff(:,1);
direction = direction / norm(direction);
normal = [-direction(2); direction(1)];

% Generate PCA line
t = linspace(min(score(:,1)), max(score(:,1)), 200);
x_pca = mu(1) + t .* direction(1);
y_pca = mu(2) + t .* direction(2);

% establish boundary for flight path
half_width = S_range/2;
x_left  = x_pca + half_width*normal(1);
y_left  = y_pca + half_width*normal(2);
x_right = x_pca - half_width*normal(1);
y_right = y_pca - half_width*normal(2);

% determine flight length of the line
line_length = max(score(:,1)) - min(score(:,1));

% plotting singular line
figure;
hold on;
measuresps('gl', 'k');

% Plot ROI overlay
shade_layer = ones(size(sampling_zone));

h_roi = imagesc(x_grid, y_grid, shade_layer);
set(gca, 'YDir', 'normal');

h_roi.AlphaData = 0.35 * double(sampling_zone);
h_roi.AlphaDataMapping = 'none';
h_roi.HandleVisibility = 'off';

colormap(gca, [1 0 0]);

% Plot PCA flight path
plot(x_pca, y_pca, 'b-', ...
     'LineWidth', 3, ...
     'HandleVisibility', 'off');

h_pca = patch([x_left fliplr(x_right)], ...
      [y_left fliplr(y_right)], ...
      'cyan', ...
      'FaceAlpha',0.3, ...
      'EdgeColor','none', ...
      'LineWidth',2);


axis image;
xlim(x_limits);
ylim(y_limits);

xlabel('x (m)');
ylabel('y (m)');

% Create proxy objects solely for the legend
h_roi_legend = patch(NaN, NaN, 'r', ...
    'FaceAlpha', 0.35, ...
    'EdgeColor', 'none');

h_gl_legend = plot(NaN, NaN, 'k-', ...
    'LineWidth', 1.5);

h_pca_legend = plot(NaN, NaN, 'b-', ...
    'LineWidth', 3);

legend([h_roi_legend, h_gl_legend, h_pca_legend], ...
    {'ROI', 'Grounding Line', 'SAR Flight Coverage'}, ...
    'Location', 'northeast');

% figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/roi_tiled_plots';
% exportgraphics(gcf, fullfile(figure_dir,'Thwaites_PCA_singular.jpg'), ...
% 'Resolution',300);

% plotting parallel lines
figure;
hold on;
measuresps('gl', 'k');

% Plot ROI overlay
shade_layer = ones(size(sampling_zone));

h_roi = imagesc(x_grid, y_grid, shade_layer);
set(gca, 'YDir', 'normal');

h_roi.AlphaData = 0.5 * double(sampling_zone);
h_roi.AlphaDataMapping = 'none';
h_roi.HandleVisibility = 'off';

colormap(gca, [1 0 0]);

plot(x_pca, y_pca, 'b-', ...
     'LineWidth', 3, ...
     'HandleVisibility', 'off');

h_pca = patch([x_left fliplr(x_right)], ...
      [y_left fliplr(y_right)], ...
      'c', ...
      'FaceAlpha',0.3, ...
      'EdgeColor','none', ...
      'LineWidth',2);

total_flight = line_length; 

% Plot PCA flight path
for i=1:1
    x_offset = x_pca + S_range*i*normal(1);
    y_offset = y_pca + S_range*i*normal(2);
    
    x_left  = x_offset + half_width*normal(1);
    y_left  = y_offset + half_width*normal(2);
    x_right = x_offset - half_width*normal(1);
    y_right = y_offset - half_width*normal(2);

    plot(x_offset, y_offset, 'b--', 'LineWidth', 2);
    patch([x_left fliplr(x_right)], ...
      [y_left fliplr(y_right)], ...
      'c', ...
      'FaceAlpha',0.3, ...
      'EdgeColor','none', ...
      'LineWidth',2);

    x_offset = x_pca - S_range*i*normal(1);
    y_offset = y_pca - S_range*i*normal(2);
    
    x_left  = x_offset + half_width*normal(1);
    y_left  = y_offset + half_width*normal(2);
    x_right = x_offset - half_width*normal(1);
    y_right = y_offset - half_width*normal(2);
    plot(x_offset, y_offset, 'b--', 'LineWidth', 2);
    
    patch([x_left fliplr(x_right)], ...
        [y_left fliplr(y_right)], ...
        'c', ...
        'FaceAlpha',0.3, ...
        'EdgeColor','none', ...
        'LineWidth',2);

    % total flight summed
    total_flight = 2*line_length*i + 2*i*S_range;
end


axis image;
xlim(x_limits);
ylim(y_limits);

xlabel('x (m)');
ylabel('y (m)');

% Create proxy objects solely for the legend
h_roi_legend = patch(NaN, NaN, 'r', ...
    'FaceAlpha', 0.35, ...
    'EdgeColor', 'none');

h_gl_legend = plot(NaN, NaN, 'k-', ...
    'LineWidth', 1.5);

h_pca_legend = plot(NaN, NaN, 'b-', ...
    'LineWidth', 3);

legend([h_roi_legend, h_gl_legend, h_pca_legend], ...
    {'ROI', 'Grounding Line', 'SAR Flight Coverage'}, ...
    'Location', 'northeast');

% figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/roi_tiled_plots';
% exportgraphics(gcf, fullfile(figure_dir,'Thwaites_PCA_parallel.jpg'), ...
% 'Resolution',300);

%% Flight time and distance analysis
% drone velocity
speed = 20*(60*60); % m/hour

% determine total flight length
fprintf('Total flight distance for area coverage = %.4f km \n', total_flight/1000);
fprintf('Total flight time for area coverage = %.4f hours \n', total_flight/speed);

% determine area of coverage of specified zone
area = length(x_gl)*grid_resolution^2; % m^2

% idealistic calculation of time required for coverage
ideal_flight_length = area / (S_range) / 1000; % km
ideal_flight_time = ideal_flight_length/(speed/1000); % hours

fprintf('Ideal distance for area coverage = %.4f km \n', ideal_flight_length);
fprintf('Ideal flight time for area coverage = %.4f hours \n', ideal_flight_time);

%% Plotting all ice streams in the Thwaites region
glacier_names = { ...
    'Pine Island', ...
    'Thwaites', ...
    'Haynes', ...
    'Pope', ...
    'Smith', ...
    'Kohler'};

colors = lines(numel(glacier_names));    % or choose your own colors

ideal_flight_lengths = cell(numel(glacier_names), 1);
sampling_zones = cell(numel(glacier_names), 1);

distance_line_handles = gobjects(0);
distance_text_handles = gobjects(0);
min_distance = NaN(numel(glacier_names));

figure;
mapzoomps('Thwaites Glacier');
hold on;
measuresps('gl', 'k', 'DisplayName', 'Grounding Line')

for i = 1:numel(glacier_names)
    % Extract grounding line for this glacier

    [xbasin,ybasin] = basin_data( ...
        'imbie refined', glacier_names{i}, 'xy');

    inside = inpolygon(xgl_all, ygl_all, xbasin, ybasin);

    xgl_glacier = xgl_all;
    ygl_glacier = ygl_all;

    xgl_glacier(~inside) = NaN;
    ygl_glacier(~inside) = NaN;

    valid = isfinite(xgl_glacier);

    gl_points = [ ...
        xgl_glacier(valid), ...
        ygl_glacier(valid)];
    
    [~,distance_to_gl] = knnsearch(gl_points,grid_points);

    sampling_zone = reshape( ...
        distance_to_gl <= radius,...
        size(X_grid));
    
    sampling_zones{i} = sampling_zone;

    % figure out idealized flight path length
    grid_size = X_grid(sampling_zone);
    area = length(grid_size)*grid_resolution^2; % m^2
    ideal_flight_lengths{i} = area / (S_range) / 1000; % km
    ideal_flight_times{i} = ideal_flight_lengths{i}/(speed/1000); % hours

    % Calculate distance between consecutive ice streams
    mask_i = sampling_zones{i};
    
    points_i = [ ...
        X_grid(mask_i), ...
        Y_grid(mask_i)];
    
    if i > 1
    
        j = i - 1;
    
        mask_j = sampling_zones{j};
    
        points_j = [ ...
            X_grid(mask_j), ...
            Y_grid(mask_j)];
    
        % For each point in region i, find nearest point in region j
        [nearest_idx_j, distances] = knnsearch(points_j, points_i);
    
        % Find the minimum of those nearest-neighbor distances
        [min_distance(i,j), idx_i] = min(distances);
    
        % Corresponding closest point in region j
        idx_j = nearest_idx_j(idx_i);
    
        % Coordinates of the two closest points
        closest_point_i = points_i(idx_i,:);
        closest_point_j = points_j(idx_j,:);
    
        % Symmetric distance matrix
        % min_distance(j,i) = min_distance(i,j);

        if min_distance(i,j) > 0
        
            min_distance_km = min_distance(i,j) / 1000;
        
            % Plot line
            h_line = plot( ...
                 [closest_point_i(1), closest_point_j(1)], ...
                 [closest_point_i(2), closest_point_j(2)], ...
                 'r-', ...
                 'LineWidth', 2, ...
                 'HandleVisibility', 'off');
            
            distance_line_handles(end+1) = h_line;
        
            % Midpoint
            x_mid = mean([closest_point_i(1), closest_point_j(1)]);
            y_mid = mean([closest_point_i(2), closest_point_j(2)]);
        
            % Unit vector along the line
            dx = closest_point_j(1) - closest_point_i(1);
            dy = closest_point_j(2) - closest_point_i(2);
            L = hypot(dx, dy);
        
            % Perpendicular unit vector
            nx = -dy / L;
            ny =  dx / L;
        
            % Offset (10 km)
            offset = 75e3;
        
            h_text = text( ...
                x_mid + offset*nx, ...
                y_mid + offset*ny, ...
                sprintf('%.0f km, %.0f minutes', min_distance_km, min_distance_km*1000/speed*60), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', ...
                'BackgroundColor', 'w', ...
                'Margin', 2, ...
                'Clipping', 'on');
            
            distance_text_handles(end+1) = h_text;
        end
    end 

    % Plot colored sampling zone
    shade = ones(size(sampling_zone));

    h = imagesc(x_grid,y_grid,shade);

    h.AlphaData = 0.35*double(sampling_zone);
    h.AlphaDataMapping = 'none';

    % Make this overlay a single solid color
    h.CData = repmat(reshape(colors(i,:),1,1,3),...
                     size(sampling_zone,1),...
                     size(sampling_zone,2));
    
    % Plot grounding line
    plot(xgl_glacier,ygl_glacier,...
        'Color',colors(i,:),...
        'LineWidth',3,...
        'DisplayName',sprintf('%s: %.0f km, %.2f hours', glacier_names{i}, ideal_flight_lengths{i}, ideal_flight_times{i}));
end

set(gca,'YDir','normal');
    axis image;
xlim(x_limits);
ylim(y_limits);

xlabel('x (m)');
ylabel('y (m)');

uistack(distance_line_handles, 'top');
uistack(distance_text_handles, 'top');
legend('Location','northeast');


% True north arrow (bottom-right)

% Arrow location (10% from right, 10% from bottom)
x0 = x_limits(2) - 0.10 * diff(x_limits);
y0 = y_limits(1) + 0.20 * diff(y_limits);

% Convert to lat/lon
[lat0, lon0] = ps2ll(x0, y0);

% Point slightly farther north
lat_north = lat0 + 0.1;
[x_north, y_north] = ll2ps(lat_north, lon0);

% Direction toward true north
dx = x_north - x0;
dy = y_north - y0;

% Normalize
L = hypot(dx,dy);
dx = dx/L;
dy = dy/L;

% Arrow length
arrow_length = 80e3;   % 80 km

dx = dx * arrow_length;
dy = dy * arrow_length;

% Draw arrow
quiver( ...
    x0, y0, ...
    dx, dy, ...
    0, ...
    'k', ...
    'LineWidth', 2.5, ...
    'MaxHeadSize', 0.8, ...
    'HandleVisibility', 'off');

% "N" label
text( ...
    x0 + 1.15*dx, ...
    y0 + 1.15*dy, ...
    'N', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 16, ...
    'FontWeight', 'bold', ...
    'BackgroundColor', 'w', ...
    'Margin', 1, ...
    'HandleVisibility', 'off');

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/roi_tiled_plots';
exportgraphics(gcf, fullfile(figure_dir,'Thwaites_multi_glaciers_tiled.jpg'), ...
'Resolution',300);