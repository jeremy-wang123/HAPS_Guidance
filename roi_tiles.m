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

% iterate through each gridpoint, determine if it switches status any 
% point within 10km radius

sampling_zone = false(size(X_grid));

for idx = 1:numel(X_grid)

    x0 = X_grid(idx);
    y0 = Y_grid(idx);

    sampling_zone(idx) = is_near_grounding_line( ...
        x0, y0, radius, search_spacing);
end

%% Grounding line plot
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

% figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/roi_tiled_plots';
% exportgraphics(gcf, fullfile(figure_dir,sprintf('Pine_Island_%dkm_tiled.jpg', radius/1000)), ...
% 'Resolution',300);

%% PCA for singular flight path

% establish the swath range
beta_el = 40; % elevation beam width
beta_az = 40; % azimuth beam width
theta = 40; % look angle
z = 20e3; % altitude of flight

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

%% plotting singular line
figure;

% Plot ROI overlay
shade_layer = ones(size(sampling_zone));

h_roi = imagesc(x_grid, y_grid, shade_layer);
set(gca, 'YDir', 'normal');

h_roi.AlphaData = 0.35 * double(sampling_zone);
h_roi.AlphaDataMapping = 'none';
h_roi.HandleVisibility = 'off';

colormap(gca, [1 0 0]);

hold on;

% Plot grounding line
measuresps('gl', 'k', 'LineWidth', 1.5);

% Plot PCA flight path
plot(x_pca, y_pca, 'b-', ...
     'LineWidth', 3, ...
     'HandleVisibility', 'off');

h_pca = patch([x_left fliplr(x_right)], ...
      [y_left fliplr(y_right)], ...
      'b', ...
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

%% plotting parallel lines
figure;

% Plot ROI overlay
shade_layer = ones(size(sampling_zone));

h_roi = imagesc(x_grid, y_grid, shade_layer);
set(gca, 'YDir', 'normal');

h_roi.AlphaData = 0.35 * double(sampling_zone);
h_roi.AlphaDataMapping = 'none';
h_roi.HandleVisibility = 'off';

colormap(gca, [1 0 0]);

hold on;

% Plot grounding line
measuresps('gl', 'k', 'LineWidth', 1.5);

plot(x_pca, y_pca, 'b-', ...
     'LineWidth', 3, ...
     'HandleVisibility', 'off');

h_pca = patch([x_left fliplr(x_right)], ...
      [y_left fliplr(y_right)], ...
      'b', ...
      'FaceAlpha',0.7, ...
      'EdgeColor','none', ...
      'LineWidth',2);

total_flight = line_length; 

% Plot PCA flight path
for i=1:5
    x_offset = x_pca + S_range*i*normal(1);
    y_offset = y_pca + S_range*i*normal(2);
    
    x_left  = x_offset + half_width*normal(1);
    y_left  = y_offset + half_width*normal(2);
    x_right = x_offset - half_width*normal(1);
    y_right = y_offset - half_width*normal(2);

    plot(x_offset, y_offset, 'b--', 'LineWidth', 2);
    patch([x_left fliplr(x_right)], ...
      [y_left fliplr(y_right)], ...
      'b', ...
      'FaceAlpha',0.7-0.1*i, ...
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
        'b', ...
        'FaceAlpha',0.7-0.1*i, ...
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
ideal_flight_time = ideal_flight_length/(speed/1000); % days

fprintf('Ideal distance for area coverage = %.4f km \n', ideal_flight_length);
fprintf('Ideal flight time for area coverage = %.4f hours \n', ideal_flight_time);

%% plot interpolated grounding line 

% Coordinate matrices
[X_grid, Y_grid] = meshgrid(x_grid, y_grid);

% Grounded/floating mask
grounded_mask = isgrounded(X_grid, Y_grid);

% Extract the grounded/floating boundary at logical level 0.5
C = contourc(x_grid, y_grid, double(grounded_mask), [0.5 0.5]);

segment_lengths = [];
segments = {};

col = 1;

while col < size(C,2)

    contour_level = C(1,col);
    n_points = C(2,col);

    % Coordinates for this contour segment
    x_segment = C(1, col+1:col+n_points);
    y_segment = C(2, col+1:col+n_points);

    % Distance between consecutive points
    dx = diff(x_segment);
    dy = diff(y_segment);

    % Length of this grounding-line segment
    segment_length = sum(hypot(dx,dy));

    segments{end+1} = [x_segment(:), y_segment(:)];
    segment_lengths(end+1) = segment_length;

    col = col + n_points + 1;
end

total_gl_length = sum(segment_lengths);

fprintf('Total grounding-line length = %.2f km\n', ...
    total_gl_length/1000);

figure;
hold on;

for k = 1:numel(segments)
    plot(segments{k}(:,1), segments{k}(:,2), ...
        'k-', 'LineWidth', 1.5);
end
measuresps('gl', 'r')
axis image;
xlim(x_limits);
ylim(y_limits);

xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Grounding-line length: %.1f km', ...
    total_gl_length/1000));