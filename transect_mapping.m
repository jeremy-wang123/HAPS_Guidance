%% Create Transect Map
clear;

%% Nearest Neighbor Approach

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
% calving_handles = scatter( ...
%     x_front, y_front, ...
%     10, 'r', 'filled', ...
%     'DisplayName', 'Calving front');

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

%% Perpendicular Transect Approach

%% Import bed machine
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

%% Manually create idealistic Thwaites shelf

% Identify ocean, floating ice, and grounded ice

ocean    = bm_mask == 0;
floating = bm_mask == 3;
grounded = bm_mask == 2;

% Thwaites outline
[wx,wy] = antbounds_data('Thwaites','xy');

inside_region = inpolygon( ...
    xgrid, ...
    ygrid, ...
    wx, ...
    wy);

% Initial floating shelf
shelf = ...
    inside_region & ...
    floating & ...
    isfinite(H);

%% ---------------------------------------------------------
% 1. Identify the REAL ocean
%
% Only ocean connected to the edge of the domain is considered
% the external ocean. This removes internal ocean holes from
% being labeled as calving fronts.
% ----------------------------------------------------------

CC_ocean = bwconncomp(ocean, 8);

external_ocean = false(size(ocean));

[nrows,ncols] = size(ocean);

for k = 1:CC_ocean.NumObjects

    inds = CC_ocean.PixelIdxList{k};

    [r,c] = ind2sub(size(ocean), inds);

    % Does this ocean component touch the edge of the domain?
    touches_edge = ...
        any(r == 1)     || ...
        any(r == nrows) || ...
        any(c == 1)     || ...
        any(c == ncols);

    if touches_edge
        external_ocean(inds) = true;
    end

end


%% ---------------------------------------------------------
% 2. Identify main grounded ice versus grounded islands
%
% Main Antarctic grounded ice should connect to the edge of
% the computational domain.
%
% Any isolated grounded component is interpreted as an island.
% ----------------------------------------------------------

% Classify grounded components

CC_grounded = bwconncomp(grounded, 8);

main_grounded = false(size(grounded));
open_ocean_island = false(size(grounded));

[nrows,ncols] = size(grounded);

% Cells adjacent to external ocean
external_ocean_neighbors = conv2( ...
    double(external_ocean), ...
    ones(3), ...
    'same') > 0;

for k = 1:CC_grounded.NumObjects

    inds = CC_grounded.PixelIdxList{k};

    [r,c] = ind2sub(size(grounded), inds);

    % Mainland grounded ice:
    % connected to edge of computational domain
    touches_domain_edge = ...
        any(r == 1)     || ...
        any(r == nrows) || ...
        any(c == 1)     || ...
        any(c == ncols);

    % Does this grounded component touch the true external ocean?
    touches_external_ocean = ...
        any(external_ocean_neighbors(inds));

    if touches_domain_edge

        % Main grounded continent
        main_grounded(inds) = true;

    elseif touches_external_ocean

        % Grounded feature on the ocean-facing side:
        % treat as part of the calving-front side
        open_ocean_island(inds) = true;

    end

    % Otherwise:
    % fully enclosed grounded island -> ignored entirely

end

%% ---------------------------------------------------------
% 3. Treat grounded islands like ocean
% ----------------------------------------------------------

effective_ocean = ...
    external_ocean | ...
    open_ocean_island;

%% ---------------------------------------------------------
% 4. Calving front
%
% Floating shelf cells adjacent to either:
%   - true external ocean
%   - grounded islands
% ----------------------------------------------------------

ocean_neighbors = conv2( ...
    double(effective_ocean), ...
    ones(3), ...
    'same') > 0;

calving_front = ...
    shelf & ...
    ocean_neighbors;


%% ---------------------------------------------------------
% 5. Grounding line
%
% ONLY shelf cells touching the main grounded ice.
% Grounded islands therefore no longer produce GL points.
% ----------------------------------------------------------

grounded_neighbors = conv2( ...
    double(main_grounded), ...
    ones(3), ...
    'same') > 0;

grounding_line = ...
    shelf & ...
    grounded_neighbors;


% Coordinates

x_front = xgrid(calving_front);
y_front = ygrid(calving_front);

x_gl = xgrid(grounding_line);
y_gl = ygrid(grounding_line);

% Force column vectors
x_front = x_front(:);
y_front = y_front(:);

x_gl = x_gl(:);
y_gl = y_gl(:);


%% Perpendicular PCA approach 

% ---------------------------------------------------------
% PCA of grounding line
% ----------------------------------------------------------

GL = [x_gl, y_gl];

% Center coordinates
GL_mean = mean(GL,1);
GL_centered = GL - GL_mean;

% PCA
[coeff, score, latent] = pca(GL);

% First principal component = strike direction
strike = coeff(:,1);

% Make sure it is a unit vector
strike = strike / norm(strike);

% Perpendicular directions
normal1 = [-strike(2); strike(1)];
normal2 = -normal1;

fprintf('Strike direction: [%.4f, %.4f]\n', ...
    strike(1), strike(2));

fprintf('Transect normal: [%.4f, %.4f]\n', ...
    normal1(1), normal1(2));
%% ---------------------------------------------------------
% Generate PCA-normal transects from calving front to GL
% ----------------------------------------------------------

% Region that transects are allowed to cross
%
% Allow:
%   - floating shelf
%   - internal ocean holes
%   - ignored internal grounded islands / mask holes
%
% Do NOT allow:
%   - external ocean

allowed_transect_region = ...
    inside_region & ...
    ~external_ocean;

n_transects = numel(x_front);

x_gl_match = NaN(n_transects,1);
y_gl_match = NaN(n_transects,1);
distance_to_gl = NaN(n_transects,1);

% Search parameters
ds = grid_spacing_m/4;     % search step
L_max = 100e3;             % maximum search distance
s = 0:ds:L_max;

% Grounding-line lookup mask
gl_mask = grounding_line;

for i = 1:n_transects

    xf = x_front(i);
    yf = y_front(i);

    best_distance = Inf;
    best_x = NaN;
    best_y = NaN;

    % Test both perpendicular directions
    normals = [normal1, normal2];

    for d = 1:2

        nvec = normals(:,d);

        % Coordinates along ray
        x_test = xf + s*nvec(1);
        y_test = yf + s*nvec(2);

        % -------------------------------------------------
        % First make sure ray does not cross external ocean
        % -------------------------------------------------

        allowed = interp2( ...
            xgrid, ...
            ygrid, ...
            double(allowed_transect_region), ...
            x_test, ...
            y_test, ...
            'nearest', ...
            0);

        % Find first invalid point
        first_bad = find(allowed < 0.5, 1);

        if isempty(first_bad)
            last_valid = numel(s);
        else
            last_valid = first_bad - 1;
        end

        if last_valid < 1
            continue
        end

        x_valid = x_test(1:last_valid);
        y_valid = y_test(1:last_valid);
        s_valid = s(1:last_valid);

        % -------------------------------------------------
        % Check for grounding-line intersection
        % -------------------------------------------------

        hits_gl = interp2( ...
            xgrid, ...
            ygrid, ...
            double(gl_mask), ...
            x_valid, ...
            y_valid, ...
            'nearest', ...
            0);

        hit_idx = find(hits_gl > 0.5, 1, 'first');

        if ~isempty(hit_idx)

            this_distance = s_valid(hit_idx);

            if this_distance < best_distance

                best_distance = this_distance;
                best_x = x_valid(hit_idx);
                best_y = y_valid(hit_idx);

            end

        end

    end

    % Save result
    if isfinite(best_distance)

        x_gl_match(i) = best_x;
        y_gl_match(i) = best_y;
        distance_to_gl(i) = best_distance;

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
    distance_to_gl > 0;

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
h_p25 = prctile(H_valid, 25, 1, 'Method', 'exact');
h_p75 = prctile(H_valid, 75, 1, 'Method', 'exact');

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

%% theoretical thickness profile
% Thickness profile parameters
h_front = mean_h_front; % m
h_gl    = mean_h_gl; % m

L = L_mean; 

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