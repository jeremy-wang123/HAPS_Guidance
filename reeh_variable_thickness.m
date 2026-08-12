%% Variable thickness profiles
clear;clc;

%% Parameter values
mean_amp = calculate_mean_amp('Thwaites', 'top2');

rho_w = 1025; % kg/m^3
g = 9.81; % m/s^2
nu = 0.325; % poisson's ratio
w0 = mean_amp; % m
sigma_z = 0.0054;  % m (vertical displacement uncertainty)

% viscoelastic rheology specific parameters
A = 3.5e-25; % Pa^-3 s^-1
tau_e = 0.1e6; % Pa 
mu_m = 1/(2*A*tau_e^2); % steady creep viscosity
E_m = 9.3e9; % young's modulus GPa
E_v = 10e9; % elastic modulus GPa
mu_v = 600e9; % elastic viscosity GPa
T = 12*3600;
omega = 2*pi/T;

%% Thickness profiles and transects
addpath('/Users/jeremywang/Documents/MATLAB/BedMachine');

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

bm_mask = bedmachine_interp( ...
    'mask', ...
    xgrid, ...
    ygrid, ...
    'method', ...
    'nearest');

ocean = bm_mask == 0;
floating = bm_mask == 3;

% Get Thwaites shelf outline
[wx,wy] = antbounds_data('Thwaites','xy');
inside_region = inpolygon(xgrid,ygrid,wx,wy);

% Define the Thwaites floating shelf
shelf = ...
    inside_region & ...
    floating & ...
    isfinite(H);

% Calving front
% Floating shelf cells adjacent to ocean
ocean_neighbors = conv2( ...
    double(ocean), ...
    ones(3), ...
    'same') > 0;

calving_front = ...
    shelf & ocean_neighbors;


% Grounding line
grounded = bm_mask == 2;

% Cells having at least one grounded neighbor
grounded_neighbors = conv2( ...
    double(grounded), ...
    ones(3), ...
    'same') > 0;

% Floating shelf cells adjacent to grounded ice
grounding_line = ...
    shelf & grounded_neighbors;

% Coordinates of front cells
x_front = xgrid(calving_front);
y_front = ygrid(calving_front);

% Coordinates of grounding-line cells
x_gl = xgrid(grounding_line);
y_gl = ygrid(grounding_line);

% Force column vectors
x_front = x_front(:);
y_front = y_front(:);

x_gl = x_gl(:);
y_gl = y_gl(:);

n_transects = numel(x_front);

% Storage
x_gl_match = NaN(n_transects,1);
y_gl_match = NaN(n_transects,1);
distance_to_gl = NaN(n_transects,1);
idx_gl = NaN(n_transects,1);

% How many nearest GL candidates to test initially
K = min(100, numel(x_gl));

% Find K nearest GL points for every front point
[idx_candidates, dist_candidates] = knnsearch( ...
    [x_gl, y_gl], ...
    [x_front, y_front], ...
    'K', K);


for i = 1:n_transects

    % Current front point
    xf = x_front(i);
    yf = y_front(i);

    % Test GL candidates from nearest to farther away
    for j = 1:K

        gl_idx = idx_candidates(i,j);

        xg = x_gl(gl_idx);
        yg = y_gl(gl_idx);

        % Straight-line distance
        L_candidate = hypot(xf - xg, yf - yg);

        % Number of samples along candidate transect
        %
        % Use spacing finer than the BedMachine grid so that
        % we don't accidentally jump across a narrow ocean gap.
        n_test = max( ...
            2, ...
            ceil(L_candidate/(grid_spacing_m/2)));

        % Straight line GL -> front
        x_test = linspace(xg, xf, n_test);
        y_test = linspace(yg, yf, n_test);

        % Determine whether every point is inside shelf
        inside_test = interp2( ...
            xgrid, ...
            ygrid, ...
            double(shelf), ...
            x_test, ...
            y_test, ...
            'nearest', ...
            0);

        % Accept this GL point only if the WHOLE straight line
        % stays inside the floating shelf
        if all(inside_test > 0.5)

            idx_gl(i) = gl_idx;

            x_gl_match(i) = xg;
            y_gl_match(i) = yg;

            distance_to_gl(i) = L_candidate;

            break
        end

    end

end

%% Build normalized thickness profile for each valid transect

n_profile = 1000;

% Normalized coordinate:
% xi = 0 -> grounding line
% xi = 1 -> calving front
xi = linspace(0, 1, n_profile);

% One row per transect
H_profiles = NaN(n_transects, n_profile);

for i = 1:n_transects

    % Skip transects where no valid GL match was found
    if ~isfinite(distance_to_gl(i)) || distance_to_gl(i) <= 0
        continue
    end

    % Coordinates along the straight GL-to-front transect
    x_line = ...
        x_gl_match(i) + ...
        xi .* (x_front(i) - x_gl_match(i));

    y_line = ...
        y_gl_match(i) + ...
        xi .* (y_front(i) - y_gl_match(i));

    % Sample BedMachine thickness along the transect
    H_profiles(i,:) = interp2( ...
        xgrid, ...
        ygrid, ...
        H, ...
        x_line, ...
        y_line, ...
        'linear');
end

%% Summary thickness quantities

% Mean thickness at detected boundaries
mean_h_front = mean(H(calving_front), 'omitnan');
mean_h_gl    = mean(H(grounding_line), 'omitnan');

% Normalized BedMachine thickness profile
n_profile = size(H_profiles, 2);

% Common normalized coordinate
% 0 = grounding line
% 1 = calving front
xi = linspace(0, 1, n_profile);

% Only keep valid transects
valid_transects = ...
    isfinite(distance_to_gl) & ...
    distance_to_gl > 0 & ...
    all(isfinite(H_profiles), 2);

H_valid = H_profiles(valid_transects, :);

% Mean profile
h_profile_mean = mean( ...
    H_valid, ...
    1, ...
    'omitnan');

% Median profile
h_profile_median = median( ...
    H_valid, ...
    1, ...
    'omitnan');

% 25th and 75th percentiles
h_p25 = prctile(H_valid, 25, 1);
h_p75 = prctile(H_valid, 75, 1);

% Mean shelf length
L_mean = mean( ...
    distance_to_gl(valid_transects), ...
    'omitnan');

%% Plot all valid transects

% Create full shelf thickness map
H_shelf_map = H;
H_shelf_map(~shelf) = NaN;

figure;
hold on;
box on;

% Shelf thickness background
pcolor(xgrid, ygrid, H_shelf_map);
shading flat;

cb = colorbar;
cb.Label.String = 'Ice thickness (m)';

% Plot ALL valid transects
for i = 1:n_transects

    if ~isfinite(distance_to_gl(i))
        continue
    end

    plot( ...
        [x_gl_match(i), x_front(i)], ...
        [y_gl_match(i), y_front(i)], ...
        '-', ...
        'Color', [0 0 0 0.3], ...   % black, 30% opacity
        'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
end

% Detected grounding-line points
gl_handles = scatter( ...
    x_gl, y_gl, ...
    10, 'k', 'filled', ...
    'DisplayName', 'Grounding line');

% Detected calving-front points
calving_handles = scatter( ...
    x_front, y_front, ...
    10, 'r', 'filled', ...
    'DisplayName', 'Calving front');

axis equal;

% Find coordinates containing valid shelf data
valid_shelf = isfinite(H_shelf_map);

x_valid = xgrid(valid_shelf);
y_valid = ygrid(valid_shelf);

% Add small padding
padding = 5e3;   % km
    
xlim([min(x_valid) - padding, ...
      max(x_valid) + padding]);

ylim([min(y_valid) - padding, ...
      max(y_valid) + padding]);
box on;


xlabel('Polar stereographic x (m)');
ylabel('Polar stereographic y (m)');

title('Shelf-Constrained GL-to-Calving-Front Transects');

legend([gl_handles, calving_handles], 'Location','southwest');
%% Determine settling time of all transects

parameters.rho_w    = rho_w;
parameters.g        = g;
parameters.E_m      = E_m;
parameters.E_v      = E_v;
parameters.mu_m     = mu_m;
parameters.mu_v     = mu_v;
parameters.nu       = nu;
parameters.w0       = w0;
parameters.omega    = omega;

% Calculate viscoelastic settling distance for every transect
n_transects = size(H_profiles, 1);
viscoelastic_settling = NaN(n_transects, 1);
elastic_settling = NaN(n_transects, 1);
viscous_settling = NaN(n_transects, 1);

for i = 1:n_transects

    % Length of this particular transect
    L_i = distance_to_gl(i);

    % Skip invalid transects
    if L_i <= 0 || any(~isfinite(H_profiles(i,:)))
        continue
    end

    % Physical x-coordinate for this transect
    x_i = xi * L_i;

    % Thickness profile for this transect
    h_i = H_profiles(i,:);

    % Parameters
    p = parameters;
    p.x_thickness = x_i;
    p.h_thickness = h_i;

    % Solve variable-thickness viscoelastic beam
    ve_i = viscoelastic_beam_profile( ...
        L_i, ...
        x_i, ...
        p);
    elastic_i = elastic_beam_profile( ...
        L_i, ...
        x_i, ...
        p);
    viscous_i = viscous_beam_profile( ...
        L_i, ...
        x_i, ...
        p);
    
    % Settling distance
    ve_settling_i = calculate_settling_distance( ...
        ve_i, ...
        x_i, ...
        w0, ...
        sigma_z); % km
    elastic_settling_i = calculate_settling_distance( ...
        elastic_i, ...
        x_i, ...
        w0, ...
        sigma_z); % km
    viscous_settling_i = calculate_settling_distance( ...
        viscous_i, ...
        x_i, ...
        w0, ...
        sigma_z); % km
    
    % If settling distance is not reached within the beam,
    % use the full beam length
    if ~isfinite(ve_settling_i) || ve_settling_i > L_i
        ve_settling_i = L_i/1e3; % km
    end
    if ~isfinite(elastic_settling_i) || elastic_settling_i > L_i
        elastic_settling_i = L_i/1e3; % km
    end
    if ~isfinite(viscous_settling_i) || viscous_settling_i > L_i
        viscous_settling_i = L_i/1e3; % km
    end

    viscoelastic_settling(i) = ve_settling_i;
    elastic_settling(i) = elastic_settling_i;
    viscous_settling(i) = viscous_settling_i;
end
%% Convert settling distances to map coordinates
[ve_x_settle, ve_y_settle, ve_valid] = map_settling_points( ...
    viscoelastic_settling, distance_to_gl, x_gl_match, y_gl_match, x_front, y_front);
[elastic_x_settle, elastic_y_settle, elastic_valid] = map_settling_points( ...
    elastic_settling, distance_to_gl, x_gl_match, y_gl_match, x_front, y_front);
[viscous_x_settle, viscous_y_settle, viscous_valid] = map_settling_points( ...
    viscous_settling, distance_to_gl, x_gl_match, y_gl_match, x_front, y_front);

%% Plot settling-distance locations

figure;
hold on;

% Ice shelf thickness
pcolor(xgrid/1e3, ygrid/1e3, H_shelf_map);
shading flat;

cb = colorbar;
cb.Label.String = 'Ice thickness (m)';

% Grounding line
h_gl_plot = plot( ...
    x_gl/1e3, ...
    y_gl/1e3, ...
    'k.', ...
    'DisplayName','Grounding line');

% Calving front
h_front_plot = plot( ...
    x_front/1e3, ...
    y_front/1e3, ...
    'r.', ...
    'DisplayName','Calving front');

% Plot ALL valid transects
for i = 1:n_transects
    if ~isfinite(distance_to_gl(i))
        continue
    end

    plot( ...
        [x_gl_match(i), x_front(i)]/1e3, ...
        [y_gl_match(i), y_front(i)]/1e3, ...
        'k-', ...
        'Color', [0 0 0 0.2], ...   % black, 30% opacity
        'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
end

%{
% Valid settling points
x_set = ve_x_settle(valid);
y_set = ve_x_settle(valid);

% Approximate Thwaites bad-point coordinates [x, y] in km
bad_points_km = [
    -1571.42, -433.159
    -1541.47, -462.574
    -1541.16, -462.874
    -1541.13, -462.946
];

bad_points = bad_points_km * 1e3;
bad_idx = knnsearch( ...
    [x_set(:), y_set(:)], ...
    bad_points);
bad_idx = unique(bad_idx);

% Print what is being removed
fprintf('Removing settling points:\n');

for i = 1:numel(bad_idx)
    fprintf('  index %d: (%.3f, %.3f) km\n', ...
        bad_idx(i), ...
        x_set(bad_idx(i))/1e3, ...
        y_set(bad_idx(i))/1e3);
end

% Remove bad points
x_set(bad_idx) = [];
y_set(bad_idx) = [];

n_set = numel(x_set);

% Order settling points by nearest-neighbor walking

% Start at one extreme point
[~, start_idx] = max(x_set);

order = NaN(n_set,1);
used = false(n_set,1);

order(1) = start_idx;
used(start_idx) = true;

for k = 2:n_set

    current_idx = order(k-1);

    % Distance from current point to all settling points
    dx = x_set - x_set(current_idx);
    dy = y_set - y_set(current_idx);

    dist = hypot(dx,dy);

    % Don't revisit points
    dist(used) = Inf;

    % Find nearest unused point
    [~, next_idx] = min(dist);

    order(k) = next_idx;
    used(next_idx) = true;
end

% Reorder settling points
x_set_ordered = x_set(order);
y_set_ordered = y_set(order);

% Distance between successive settling points
d_set = hypot( ...
    diff(x_set_ordered), ...
    diff(y_set_ordered));

% Break connection if settling points are too far apart
max_gap = 20e3;   %

x_plot = x_set_ordered;
y_plot = y_set_ordered;

break_idx = find(d_set > max_gap);

for k = fliplr(break_idx(:)')

    x_plot = [ ...
        x_plot(1:k); ...
        NaN; ...
        x_plot(k+1:end)];

    y_plot = [ ...
        y_plot(1:k); ...
        NaN; ...
        y_plot(k+1:end)];
end

plot( ...
    x_plot/1e3, ...
    y_plot/1e3, ...
    '-', ...
    'Color', [1 0.5 0], ...
    'LineWidth', 2, ...
    'DisplayName', 'Viscoelastic settling boundary');

%}

% Settling points
ve_settle = scatter( ...
    ve_x_settle(ve_valid)/1e3, ...
    ve_y_settle(ve_valid)/1e3, ...
    10, ...
    [1 0.5 0], ...
    'filled', ...
    'DisplayName','Viscoelastic settling distance');

% Elastic points
elastic_settle = scatter( ...
    elastic_x_settle(elastic_valid)/1e3, ...
    elastic_y_settle(elastic_valid)/1e3, ...
    10, ...
    'red', ...
    'filled', ...
    'DisplayName','Elastic settling distance');

% viscous points
viscous_settle = scatter( ...
    viscous_x_settle(viscous_valid)/1e3, ...
    viscous_y_settle(viscous_valid)/1e3, ...
    10, ...
    'green', ...
    'filled', ...
    'DisplayName','Viscous settling distance');

axis equal;

% Find coordinates containing valid shelf data
valid_shelf = isfinite(H_shelf_map);

x_valid = xgrid(valid_shelf)/1e3;
y_valid = ygrid(valid_shelf)/1e3;

% Add small padding
padding = 5;   % km

xlim([min(x_valid) - padding, ...
      max(x_valid) + padding]);

ylim([min(y_valid) - padding, ...
      max(y_valid) + padding]);
box on;

xlabel('x (km)');
ylabel('y (km)');

title('Thwaites ROI');

legend( ...
    [h_gl_plot, h_front_plot, ve_settle, elastic_settle, viscous_settle], ...
    'Location','southwest');

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/variable_thickness';
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end
exportgraphics(gcf, fullfile(figure_dir,'thwaites_roi_comparison.jpg'), ...
'Resolution',300);

%% theoretical thickness profile
% Thickness profile parameters
h_front = mean_h_front; % m
h_gl    = mean_h_gl; % m

% flexurarl length
I_analytical = h_front^3 / 12;          % second moment per unit width, m^3
D_analytical = E_m*I_analytical / (1 - nu^2);  % flexural rigidity per unit width, N m
lambda = (rho_w*g/(4*D_analytical))^(1/4); 
L = 30/lambda; 

x_thickness = linspace(0, L, 1000);
x_thickness_km = x_thickness/1e3;

% uniform profile
h_uniform =  h_front * ones(size(x_thickness));

% Linear profile
h_linear = ...
    h_gl + (h_front - h_gl)*(x_thickness/L);

% Quadratic taper
p_shape = 2;
h_quad = ...
    h_gl + ...
    (h_front - h_gl)*(x_thickness/L).^p_shape;

% Exponential taper
k = 2;  % controls curvature
h_exp = h_gl + ...
    (h_front - h_gl) .* ...
    (1 - exp(-k*x_thickness/L)) ./ ...
    (1 - exp(-k));

% x coords for normalized
x_norm = xi * L;
x_norm_km = x_norm/1e3; 

% Plot
figure;
hold on;
grid on;
box on;

fill( ...
    [x_norm_km fliplr(x_norm_km)], ...
    [h_p25 fliplr(h_p75)], ...
    [0.8 0.8 0.8], ...
    'EdgeColor','none', ...
    'DisplayName','25th–75th percentile');

plot(x_norm_km, h_profile_mean, ...
    'k', ...
    'LineWidth',2, ...
    'DisplayName','Normalized Mean');

plot(x_thickness_km, h_linear, ...
    'LineWidth', 2, ...
    'DisplayName', 'Linear');

plot(x_thickness_km, h_quad, ...
    'LineWidth', 2, ...
    'DisplayName', 'Quadratic');

plot(x_thickness_km, h_exp, ...
    'LineWidth', 2, ...
    'DisplayName', 'Exponential');

xlabel('Distance from grounding line (km)');
ylabel('Ice thickness (m)');
title('Ice Thickness Profiles');
set(gca, 'YDir', 'reverse');
legend('Location','best');

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/variable_thickness';
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end
% exportgraphics(gcf, fullfile(figure_dir,'elastic_thickness_profiles.jpg'), ...
% 'Resolution',300);

%% creating thickness profiles
x = linspace(0, L, 1000);
x_km = x/1e3;

% uniform profile
uniform_parameters.rho_w = rho_w;
uniform_parameters.g     = g;
uniform_parameters.E_m   = E_m;
uniform_parameters.E_v   = E_v;
uniform_parameters.mu_m  = mu_m;
uniform_parameters.mu_v  = mu_v;
uniform_parameters.nu    = nu;
uniform_parameters.w0    = w0;
uniform_parameters.omega    = omega;

uniform_parameters.x_thickness = x_thickness;
uniform_parameters.h_thickness = h_uniform;

% linear profile
linear_parameters.rho_w     = rho_w;
linear_parameters.g         = g;
linear_parameters.E_m      = E_m;
linear_parameters.E_v      = E_v;
linear_parameters.mu_m     = mu_m;
linear_parameters.mu_v     = mu_v;
linear_parameters.nu        = nu;
linear_parameters.w0        = w0;
linear_parameters.omega    = omega;

linear_parameters.x_thickness = x_thickness;
linear_parameters.h_thickness = h_quad;

% quad profile
quad_parameters.rho_w    = rho_w;
quad_parameters.g        = g;
quad_parameters.E_m      = E_m;
quad_parameters.E_v      = E_v;
quad_parameters.mu_m     = mu_m;
quad_parameters.mu_v     = mu_v;
quad_parameters.nu       = nu;
quad_parameters.w0       = w0;
quad_parameters.omega    = omega;

quad_parameters.x_thickness = x_thickness;
quad_parameters.h_thickness = h_linear;

% exponential profile
exp_parameters.rho_w    = rho_w;
exp_parameters.g        = g;
exp_parameters.E_m      = E_m;
exp_parameters.E_v      = E_v;
exp_parameters.mu_m     = mu_m;
exp_parameters.mu_v     = mu_v;
exp_parameters.nu       = nu;
exp_parameters.w0       = w0;
exp_parameters.omega    = omega;

exp_parameters.x_thickness = x_thickness;
exp_parameters.h_thickness = h_exp;

% Thwaites normalized mean thickness profile
normalized_parameters.rho_w    = rho_w;
normalized_parameters.g        = g;
normalized_parameters.E_m      = E_m;
normalized_parameters.E_v      = E_v;
normalized_parameters.mu_m     = mu_m;
normalized_parameters.mu_v     = mu_v;
normalized_parameters.nu       = nu;
normalized_parameters.w0       = w0;
normalized_parameters.omega    = omega;

normalized_parameters.x_thickness = x_norm;
normalized_parameters.h_thickness = h_profile_mean;

%% elastic beam thickness comparison
elastic_linear = elastic_beam_profile(L, x, linear_parameters);
elastic_quad = elastic_beam_profile(L, x, quad_parameters);
elastic_exp = elastic_beam_profile(L, x, exp_parameters);
elastic_norm = elastic_beam_profile(L, x, normalized_parameters);

% elastic analytical
w_uniform = w0 .* ...
    (1 - exp(-lambda*x) .* ...
    (cos(lambda*x) + sin(lambda*x)));

% Plot displacement for elastic
figure;
hold on;
grid on;
box on;

plot(x_km, elastic_linear, ...
    'LineWidth', 2, ...
    'DisplayName', 'Linear');

plot(x_km, elastic_quad, ...
    'LineWidth', 2, ...
    'DisplayName', 'Quadratic');

plot(x_km, elastic_exp, ...
    'LineWidth', 2, ...
    'DisplayName', 'Exponential');

plot(x_km, w_uniform, '--', ...
    'LineWidth', 2, ...
    'DisplayName', sprintf('Uniform: %.0f m', h_front));

plot(x_km, elastic_norm, ...
    'LineWidth', 2, ...
    'DisplayName', 'Normalized Mean');

% Upper and lower uncertainty bounds around steady state
upper_bound = w0 + sigma_z;
lower_bound = w0 - sigma_z;

yline(upper_bound, ':', ...
    '$\pm\sigma_z$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right', ...
    'HandleVisibility', 'off');
yline(lower_bound, ':', ...
    'HandleVisibility', 'off');

% settling time
uniform_settling = ...
    calculate_settling_distance( ...
        w_uniform, x, w0, sigma_z);

linear_settling = ...
    calculate_settling_distance( ...
        elastic_linear, x, w0, sigma_z);

quad_settling = ...
    calculate_settling_distance( ...
        elastic_quad, x, w0, sigma_z);

exp_settling = ...
    calculate_settling_distance( ...
        elastic_exp, x, w0, sigma_z);

norm_settling = ...
    calculate_settling_distance( ...
        elastic_norm, x, w0, sigma_z);

xline(uniform_settling, ':', ...
    sprintf('Uniform = %.1f km', uniform_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(linear_settling, ':', ...
    sprintf('Linear = %.1f km', linear_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(quad_settling, ':', ...
    sprintf('Quad = %.1f km', quad_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(exp_settling, ':', ...
    sprintf('Exp = %.1f km', exp_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(norm_settling, ':', ...
    sprintf('Norm = %.1f km', norm_settling), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(L_mean/1e3, '--', ...
    sprintf('Mean GL Distance = %.1f km', L_mean/1e3), ...
    'Color', 'red', ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'right', ...
    'HandleVisibility', 'off');

yline(w0, ':', ...
    '$w_0$', ...
    'LabelHorizontalAlignment','left', ...
    'Interpreter', 'latex', ...
    'HandleVisibility', 'off');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title(sprintf('Elastic uniform thickness vs non-uniform: h front = %.0f m, h gl = %.0f km', ...
    h_front, h_gl));

legend('Location', 'southeast');


figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/variable_thickness';
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end
% exportgraphics(gcf, fullfile(figure_dir,'elastic_thickness_comparison.jpg'), ...
% 'Resolution',300);

%% viscous beam thickness comparison
viscous_uniform = viscous_beam_profile(L, x, uniform_parameters);
viscous_linear = viscous_beam_profile(L, x, linear_parameters);
viscous_quad = viscous_beam_profile(L, x, quad_parameters);
viscous_exp = viscous_beam_profile(L, x, exp_parameters);
viscous_norm = viscous_beam_profile(L, x, normalized_parameters);

% Plot displacement for elastic
figure;
hold on;
grid on;
box on;

plot(x_km, viscous_linear, ...
    'LineWidth', 2, ...
    'DisplayName', 'Linear');

plot(x_km, viscous_quad, ...
    'LineWidth', 2, ...
    'DisplayName', 'Quadratic');

plot(x_km, viscous_exp, ...
    'LineWidth', 2, ...
    'DisplayName', 'Exponential');

plot(x_km, viscous_uniform, '--', ...
    'LineWidth', 2, ...
    'DisplayName', sprintf('Uniform: %.0f m', h_front));

plot(x_km, viscous_norm, ...
    'LineWidth', 2, ...
    'DisplayName', 'Normalized Mean');

% Upper and lower uncertainty bounds around steady state
upper_bound = w0 + sigma_z;
lower_bound = w0 - sigma_z;

yline(upper_bound, ':', ...
    '$\pm\sigma_z$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right', ...
    'HandleVisibility', 'off');
yline(lower_bound, ':', ...
    'HandleVisibility', 'off');

% settling time
uniform_settling = ...
    calculate_settling_distance( ...
        viscous_uniform, x, w0, sigma_z);

linear_settling = ...
    calculate_settling_distance( ...
        viscous_linear, x, w0, sigma_z);

quad_settling = ...
    calculate_settling_distance( ...
        viscous_quad, x, w0, sigma_z);

exp_settling = ...
    calculate_settling_distance( ...
        viscous_exp, x, w0, sigma_z);

norm_settling = ...
    calculate_settling_distance( ...
        viscous_norm, x, w0, sigma_z);

xline(uniform_settling, ':', ...
    sprintf('Uniform = %.1f km', uniform_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(linear_settling, ':', ...
    sprintf('Linear = %.1f km', linear_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(quad_settling, ':', ...
    sprintf('Quad = %.1f km', quad_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(exp_settling, ':', ...
    sprintf('Exp = %.1f km', exp_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(norm_settling, ':', ...
    sprintf('Norm = %.1f km', norm_settling), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(L_mean/1e3, '--', ...
    sprintf('Mean GL Distance = %.1f km', L_mean/1e3), ...
    'Color', 'red', ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'right', ...
    'HandleVisibility', 'off');

yline(w0, ':', ...
    '$w_0$', ...
    'LabelHorizontalAlignment','left', ...
    'Interpreter', 'latex', ...
    'HandleVisibility', 'off');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title(sprintf('Viscous uniform thickness vs non-uniform: h front = %.0f m, h gl = %.0f km', ...
    h_front, h_gl));

legend('Location', 'southeast');


figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/variable_thickness';
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end
% exportgraphics(gcf, fullfile(figure_dir,'viscous_thickness_comparison.jpg'), ...
% 'Resolution',300);

%% viscoelastic beam thickness comparison
ve_uniform = viscoelastic_beam_profile(L, x, uniform_parameters);
ve_linear = viscoelastic_beam_profile(L, x, linear_parameters);
ve_quad = viscoelastic_beam_profile(L, x, quad_parameters);
ve_exp = viscoelastic_beam_profile(L, x, exp_parameters);
ve_norm = viscoelastic_beam_profile(L, x, normalized_parameters);

% Plot displacement for elastic
figure;
hold on;
grid on;
box on;

plot(x_km, ve_linear, ...
    'LineWidth', 2, ...
    'DisplayName', 'Linear');

plot(x_km, ve_quad, ...
    'LineWidth', 2, ...
    'DisplayName', 'Quadratic');

plot(x_km, ve_exp, ...
    'LineWidth', 2, ...
    'DisplayName', 'Exponential');

plot(x_km, ve_uniform, '--', ...
    'LineWidth', 2, ...
    'DisplayName', sprintf('Uniform: %.0f m', h_front));

plot(x_km, ve_norm, ...
    'LineWidth', 2, ...
    'DisplayName', 'Normalized Mean');

% Upper and lower uncertainty bounds around steady state
upper_bound = w0 + sigma_z;
lower_bound = w0 - sigma_z;

yline(upper_bound, ':', ...
    '$\pm\sigma_z$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right', ...
    'HandleVisibility', 'off');
yline(lower_bound, ':', ...
    'HandleVisibility', 'off');

% settling time
uniform_settling = ...
    calculate_settling_distance( ...
        ve_uniform, x, w0, sigma_z);

linear_settling = ...
    calculate_settling_distance( ...
        ve_linear, x, w0, sigma_z);

quad_settling = ...
    calculate_settling_distance( ...
        ve_quad, x, w0, sigma_z);

exp_settling = ...
    calculate_settling_distance( ...
        ve_exp, x, w0, sigma_z);

norm_settling = ...
    calculate_settling_distance( ...
        ve_norm, x, w0, sigma_z);

xline(uniform_settling, ':', ...
    sprintf('Uniform = %.1f km', uniform_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(linear_settling, ':', ...
    sprintf('Linear = %.1f km', linear_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(quad_settling, ':', ...
    sprintf('Quad = %.1f km', quad_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(exp_settling, ':', ...
    sprintf('Exp = %.1f km', exp_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(norm_settling, ':', ...
    sprintf('Norm = %.1f km', norm_settling), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(L_mean/1e3, '--', ...
    sprintf('Mean GL Distance = %.1f km', L_mean/1e3), ...
    'Color', 'red', ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'right', ...
    'HandleVisibility', 'off');

yline(w0, ':', ...
    '$w_0$', ...
    'LabelHorizontalAlignment','left', ...
    'Interpreter', 'latex', ...
    'HandleVisibility', 'off');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title(sprintf('Viscoelastic uniform thickness vs non-uniform: h front = %.0f m, h gl = %.0f km', ...
    h_front, h_gl));

legend('Location', 'southeast');


figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/variable_thickness';
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end
% exportgraphics(gcf, fullfile(figure_dir,'viscoelastic_thickness_comparison.jpg'), ...
% 'Resolution',300);

%% Comparisons between elastic, viscous, and viscoelastic
% Plot displacement for elastic
figure;
hold on;
grid on;
box on;

plot(x_km, elastic_norm, ...
    'LineWidth', 2, ...
    'DisplayName', 'elastic');

plot(x_km, viscous_norm, ...
    'LineWidth', 2, ...
    'DisplayName', 'viscous');

plot(x_km, ve_norm, ...
    'LineWidth', 2, ...
    'DisplayName', 'viscoelastic');

% Upper and lower uncertainty bounds around steady state
upper_bound = w0 + sigma_z;
lower_bound = w0 - sigma_z;

yline(upper_bound, ':', ...
    '$\pm\sigma_z$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right', ...
    'HandleVisibility', 'off');
yline(lower_bound, ':', ...
    'HandleVisibility', 'off');

% settling time
elastic_settling = ...
    calculate_settling_distance( ...
        elastic_norm, x, w0, sigma_z);

viscous_settling = ...
    calculate_settling_distance( ...
        viscous_norm, x, w0, sigma_z);

ve_settling = ...
    calculate_settling_distance( ...
        ve_norm, x, w0, sigma_z);

xline(elastic_settling, ':', ...
    sprintf('Elastic = %.1f km', elastic_settling), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(viscous_settling, ':', ...
    sprintf('Viscous = %.1f km', viscous_settling), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(ve_settling, ':', ...
    sprintf('Viscoelastic = %.1f km', ve_settling), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(L_mean/1e3, '--', ...
    sprintf('Mean GL Distance = %.1f km', L_mean/1e3), ...
    'Color', 'red', ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'right', ...
    'HandleVisibility', 'off');

yline(w0, ':', ...
    '$w_0$', ...
    'LabelHorizontalAlignment','left', ...
    'Interpreter', 'latex', ...
    'HandleVisibility', 'off');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title(sprintf('Comparison uniform thickness vs non-uniform: h front = %.0f m, h gl = %.0f km', ...
    h_front, h_gl));

legend('Location', 'southeast');


figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/variable_thickness';
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end
% exportgraphics(gcf, fullfile(figure_dir,'rheology_comparison.jpg'), ...
% 'Resolution',300);

%% Local functions 
function [x_settle, y_settle, valid] = map_settling_points( ...
    settling_distance, distance_to_gl, ...
    x_gl, y_gl, x_front, y_front)

% MAP_SETTLING_POINTS
% Converts settling distances along transects into map coordinates.
%
% Inputs:
%   settling_distance - settling distance along each transect [km]
%   distance_to_gl    - distance from grounding line to front [m]
%   x_gl, y_gl        - grounding-line coordinates [m]
%   x_front, y_front  - calving-front coordinates [m]
%
% Outputs:
%   x_settle, y_settle - settling-point map coordinates [m]
%   valid              - logical array of valid settling points

    % Make everything column vectors
    settling_distance = settling_distance(:);
    distance_to_gl    = distance_to_gl(:);
    x_gl              = x_gl(:);
    y_gl              = y_gl(:);
    x_front           = x_front(:);
    y_front           = y_front(:);

    % Fraction of distance from grounding line to front
    settling_fraction = ...
        settling_distance * 1e3 ./ distance_to_gl;

    % Only keep settling points that fall along the transect
    valid = ...
        isfinite(settling_fraction) & ...
        settling_fraction >= 0 & ...
        settling_fraction <= 1;

    % Initialize outputs
    x_settle = NaN(size(settling_fraction));
    y_settle = NaN(size(settling_fraction));

    % Interpolate from grounding line toward calving front
    x_settle(valid) = ...
        x_gl(valid) + ...
        settling_fraction(valid) .* ...
        (x_front(valid) - x_gl(valid));

    y_settle(valid) = ...
        y_gl(valid) + ...
        settling_fraction(valid) .* ...
        (y_front(valid) - y_gl(valid));

end

%% Local functions for elastic beam

function w_elastic = elastic_beam_profile(L, x, parameters)
%%outputs elastic profile for varying bed thickness
% Inputs:
% L (length of beam)
% x (horizontal evaluation points)
% parameters (includes rho_w, g, E_m, nu, w0, x_thickness, h_thickness)
% Outputs:
% w_elastic (displacement profile)

    % Initial mesh and initial guess
    x_mesh = linspace(0, L, 300);

    solinit = bvpinit( ...
    x_mesh, ...
    @(x) initial_guess(x, L, parameters.w0));

    % Solve boundary-value problem
    options = bvpset( ...
        'RelTol', 1e-7, ...
        'AbsTol', 1e-9, ...
        'Stats', 'on');

    sol = bvp4c( ...
    @(x,y) beam_ode(x, y, parameters), ...
    @beam_boundary_conditions, ...
    solinit, ...
    options);

    Y = deval(sol, x);

    w_elastic     = Y(1,:);
    theta = Y(2,:);
    M     = Y(3,:);
    Q     = Y(4,:);

end

function dydx = beam_ode(x, y, p)
    w     = y(1);
    theta = y(2);
    M     = y(3);
    Q     = y(4);

    % Local thickness
    h_local = interp1( ...
        p.x_thickness, ...
        p.h_thickness, ...
        x, ...
        'linear');

    % Local second moment of area
    I_local = h_local^3/12;

    % Local flexural rigidity
    D_local = ...
        p.E_m*I_local/(1 - p.nu^2);

    dydx = [
        theta
        M/D_local
        Q
        p.rho_w*p.g*(p.w0 - w)
    ];
end

function residual = beam_boundary_conditions(ya, yb)
% At grounding line, x = 0:
% w = 0
% theta = 0
%
% At free terminus, x = L:
% M = 0
% Q = 0
    residual = [
        ya(1)
        ya(2)
        yb(3)
        yb(4)
    ];
end


function y_guess = initial_guess(x, L, w0)
% Smooth initial guess for bvp4c

    transition_length = 0.2*L;

    w_guess = w0*(1 - exp(-x/transition_length));

    theta_guess = ...
        (w0/transition_length)*exp(-x/transition_length);

    y_guess = [
        w_guess
        theta_guess
        0
        0
    ];
end

%% Local functions for viscous beam
function u_viscous_peak = viscous_beam_profile(L, x, parameters)
%%outputs viscous profile for varying bed thickness
% Inputs:
% L (length of beam)
% x (horizontal evaluation points)
% parameters (includes rho_w, g, E_m, nu, w0, x_thickness, h_thickness)
% Outputs:
% u_viscous_peak (displacement profile)
    % Initial guess
    x_mesh = linspace(0, L, 300);
    
    solinit = bvpinit( ...
        x_mesh, ...
        @(x) viscous_initial_guess(x, L, parameters.w0));
    
    % Solve
    options = bvpset( ...
        'RelTol', 1e-7, ...
        'AbsTol', 1e-9, ...
        'NMax', 20000, ...
        'Stats', 'on');
    
    sol = bvp4c( ...
        @(x,y) viscous_beam_ode(x, y, parameters), ...
        @viscous_beam_bc, ...
        solinit, ...
        options);
    
    Y = deval(sol, x);
    x_km = x/1e3;
    
    us      = Y(1,:);
    uc      = Y(2,:);
    alpha_s = Y(3,:);
    alpha_c = Y(4,:);
    Qs      = Y(5,:);
    Qc      = Y(6,:);
    Ms      = Y(7,:);
    Mc      = Y(8,:);
    
    displacement_amplitude = hypot(us, uc);
    tilt_amplitude = hypot(alpha_s, alpha_c);
    
    displacement_phase = atan2d(uc, us);
    tilt_phase = atan2d(alpha_c, alpha_s);
    
    t_peak = pi/(2*parameters.omega);
    
    u_viscous_peak = ...
        us .* sin(parameters.omega*t_peak) + ...
        uc .* cos(parameters.omega*t_peak);
end

function dydx = viscous_beam_ode(x, y, p)
% Reeh et al. four-element viscoelastic beam.
%
% State vector:
% y(1) = us       sine displacement component
% y(2) = uc       cosine displacement component
% y(3) = alpha_s  sine tilt component
% y(4) = alpha_c  cosine tilt component
% y(5) = Qs       sine shear component
% y(6) = Qc       cosine shear component
% y(7) = Ms       sine moment component
% y(8) = Mc       cosine moment component

    us      = y(1);
    uc      = y(2);
    alpha_s = y(3);
    alpha_c = y(4);
    Qs      = y(5);
    Qc      = y(6);
    Ms      = y(7);
    Mc      = y(8);

    % Local thickness
    h_local = interp1( ...
        p.x_thickness, ...
        p.h_thickness, ...
        x, ...
        'linear', ...
        'extrap');

    % Local second moment of area
    I_local = h_local^3/12;

    % Material coefficients
    omega = p.omega;
    mu_m = p.mu_m; 
    Dv_local = 4*I_local*mu_m;


    % Return derivatives of all eight state variables
    dydx = [
        alpha_s                         % dus/dx
        alpha_c                         % duc/dx
        Mc/(Dv_local*omega)             % dalpha_s/dx
       -Ms/(Dv_local*omega)             % dalpha_c/dx
        p.rho_w*p.g*(p.w0 - us)         % dQs/dx
       -p.rho_w*p.g*uc                  % dQc/dx
        Qs                              % dMs/dx
        Qc                              % dMc/dx
    ];
end

function residual = viscous_beam_bc(ya, yb)
    residual = [
        ya(1)    % us(0) = 0
        ya(2)    % uc(0) = 0
        ya(3)    % alpha_s(0) = 0
        ya(4)    % alpha_c(0) = 0
        yb(5)    % Qs(L) = 0
        yb(6)    % Qc(L) = 0
        yb(7)    % Ms(L) = 0
        yb(8)    % Mc(L) = 0
    ];
end

function y_guess = viscous_initial_guess(x, L, a)

    transition_length = 0.20*L;

    us_guess = a*(1 - exp(-x/transition_length));

    alpha_s_guess = ...
        (a/transition_length)*exp(-x/transition_length);

    y_guess = [
        us_guess
        0
        alpha_s_guess
        0
        0
        0
        0
        0
    ];
end

%% Local functions for viscoelastic beam 
function u_viscoelastic_peak = viscoelastic_beam_profile(L, x, parameters)
%%outputs u_viscoelastic_peak profile for varying bed thickness
% Inputs:
% L (length of beam)
% x (horizontal evaluation points)
% parameters (includes rho_w, g, E_m, nu, w0, x_thickness, h_thickness)
% Outputs:
% u_viscoelastic_peak (displacement profile)
    % Initial guess
    x_mesh = linspace(0, L, 300);
    
    solinit = bvpinit( ...
        x_mesh, ...
        @(x) viscoelastic_initial_guess(x, L, parameters.w0));
    
    % Solve
    options = bvpset( ...
        'RelTol', 1e-7, ...
        'AbsTol', 1e-9, ...
        'NMax', 20000, ...
        'Stats', 'on');
    
    sol = bvp4c( ...
        @(x,y) viscoelastic_beam_ode(x, y, parameters), ...
        @viscoelastic_beam_bc, ...
        solinit, ...
        options);

    Y = deval(sol, x);
    x_km = x/1e3;
    
    us      = Y(1,:);
    uc      = Y(2,:);
    alpha_s = Y(3,:);
    alpha_c = Y(4,:);
    Qs      = Y(5,:);
    Qc      = Y(6,:);
    Ms      = Y(7,:);
    Mc      = Y(8,:);
    
    displacement_amplitude = hypot(us, uc);
    tilt_amplitude = hypot(alpha_s, alpha_c);
    
    displacement_phase = atan2d(uc, us);
    tilt_phase = atan2d(alpha_c, alpha_s);
    
    t_peak = pi/(2*parameters.omega);
    
    u_viscoelastic_peak = ...
        us .* sin(parameters.omega*t_peak) + ...
        uc .* cos(parameters.omega*t_peak);
end 

function dydx = viscoelastic_beam_ode(x, y, p)
% Reeh et al. four-element viscoelastic beam.
%
% State vector:
% y(1) = us       sine displacement component
% y(2) = uc       cosine displacement component
% y(3) = alpha_s  sine tilt component
% y(4) = alpha_c  cosine tilt component
% y(5) = Qs       sine shear component
% y(6) = Qc       cosine shear component
% y(7) = Ms       sine moment component
% y(8) = Mc       cosine moment component

    us      = y(1);
    uc      = y(2);
    alpha_s = y(3);
    alpha_c = y(4);
    Qs      = y(5);
    Qc      = y(6);
    Ms      = y(7);
    Mc      = y(8);
    
    % Local thickness
    h_local = interp1( ...
        p.x_thickness, ...
        p.h_thickness, ...
        x, ...
        'linear', ...
        'extrap');

    % Local second moment of area
    I_local = h_local^3/12;

    % Material coefficients
    C1 = 3*p.mu_v/p.E_v;

    C2 = 2.25*p.mu_v / ...
        (p.E_v*p.E_m*I_local);

    C3 = 3/(4*I_local) * ...
        (1/p.E_v + ...
         1/p.E_m + ...
         p.mu_v/(p.mu_m*p.E_v));

    C4 = 2.25/(p.mu_m*I_local);

    omega = p.omega;

    % Solve the two coupled curvature equations

    curvature_matrix = [
        omega^2*C1, omega
        omega,      omega^2*C1
    ];

    curvature_rhs = [
        (omega^2*C2 + C4)*Ms + omega*C3*Mc
        omega*C3*Ms + (omega^2*C2 + C4)*Mc
    ];

    curvature = curvature_matrix \ curvature_rhs;

    dalpha_s_dx = curvature(1);
    dalpha_c_dx = curvature(2);

    % Return derivatives of all eight state variables
    dydx = [
        alpha_s                         % dus/dx
        alpha_c                         % duc/dx
        dalpha_s_dx                     % dalpha_s/dx
        dalpha_c_dx                     % dalpha_c/dx
        p.rho_w*p.g*(p.w0 - us)          % dQs/dx
       -p.rho_w*p.g*uc                  % dQc/dx
        Qs                              % dMs/dx
        Qc                              % dMc/dx
    ];
end

function residual = viscoelastic_beam_bc(ya, yb)

    residual = [
        ya(1)    % us(0) = 0
        ya(2)    % uc(0) = 0
        ya(3)    % alpha_s(0) = 0
        ya(4)    % alpha_c(0) = 0
        yb(5)    % Qs(L) = 0
        yb(6)    % Qc(L) = 0
        yb(7)    % Ms(L) = 0
        yb(8)    % Mc(L) = 0
    ];
end

function y_guess = viscoelastic_initial_guess(x, L, a)

    transition_length = 0.20*L;

    us_guess = a*(1 - exp(-x/transition_length));

    alpha_s_guess = ...
        (a/transition_length)*exp(-x/transition_length);

    y_guess = [
        us_guess
        0
        alpha_s_guess
        0
        0
        0
        0
        0
    ];
end