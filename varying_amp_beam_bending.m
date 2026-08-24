%% Varying Amplitude Beam Bending
clear;clc;
figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/variable_amplitude';
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

%% Parameter values

rho_w = 1025; % kg/m^3
g = 9.81; % m/s^2
nu = 0.325; % poisson's ratio
sigma_z = 0.0054;  % m (vertical displacement uncertainty)

% viscoelastic rheology specific parameters
A = 3.5e-25; % Pa^-3 s^-1
tau_e = 0.1e6; % Pa 
mu_m = 1/(2*A*tau_e^2); % steady creep viscosity
E_m = 9.3e9; % young's modulus GPa
E_v = 10e9; % elastic modulus GPa
mu_v = 600e9; % elastic viscosity GPa
T = 24*3600;
omega = 2*pi/T;

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

% ---------------------------------------------------------
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

% ---------------------------------------------------------
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

% ---------------------------------------------------------
% 3. Treat grounded islands like ocean
% ----------------------------------------------------------

effective_ocean = ...
    external_ocean | ...
    open_ocean_island;

% ---------------------------------------------------------
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


% ---------------------------------------------------------
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
    H_i = interp2( ...
        xgrid, ...
        ygrid, ...
        H, ...
        x_line, ...
        y_line, ...
        'linear');

    % Interpolate across internal NaN gaps
    H_i = fillmissing(H_i, 'linear');

    % Store completed thickness profile
    H_profiles(i,:) = H_i;

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

%% varying amplitude

% creating thickness profiles
x = linspace(0, L, 1000);
x_km = x/1e3;

% Thwaites normalized mean thickness profile
parameters.rho_w    = rho_w;
parameters.g        = g;
parameters.E_m      = E_m;
parameters.E_v      = E_v;
parameters.mu_m     = mu_m;
parameters.mu_v     = mu_v;
parameters.nu       = nu;
parameters.omega    = omega;

parameters.x_thickness = x_norm;
parameters.h_thickness = h_profile_mean;

% solve profiles
model_names = {'Elastic', 'Viscous', 'Maxwell'};

for j = 1:3

    figure;
    hold on;
    grid on;
    box on;

    for i = 1:10

        % Mean amplitude using top i tidal constituents
        parameters.w0 = calculate_mean_amp('Thwaites', i);

        % Calculate profile depending on model
        switch j
            case 1
                w = elastic_beam_profile(L, x, parameters);

            case 2
                w = viscous_beam_profile(L, x, parameters);

            case 3
                w = maxwell_beam_profile(L, x, parameters);
        end

        % Plot displacement
        plot(x_km, w, ...
            'LineWidth', 2, ...
            'DisplayName', sprintf('Top %d', i));

        settling_distance = ...
            calculate_settling_distance( ...
                w, x, parameters.w0, sigma_z);
        
        xline(settling_distance, ':', ...
            sprintf('Top %d = %.1f km', i, settling_distance), ...
            'LabelVerticalAlignment', 'bottom', ...
            'LabelHorizontalAlignment', 'left', ...
            'HandleVisibility', 'off');
    end

    % Mean GL distance
    xline(L_mean/1e3, '--', ...
        sprintf('Mean GL Distance = %.1f km', L_mean/1e3), ...
        'Color', 'red', ...
        'LabelVerticalAlignment', 'middle', ...
        'LabelHorizontalAlignment', 'right', ...
        'HandleVisibility', 'off');

    xlabel('Distance from grounding line (km)');
    ylabel('Vertical displacement (m)');

    title(sprintf('%s Beam Response', model_names{j}));

    legend('Location', 'southeast');

end

%%

% creating thickness profiles
x = linspace(0, L, 1000);
x_km = x/1e3;

% Thwaites normalized mean thickness profile
parameters.rho_w    = rho_w;
parameters.g        = g;
parameters.E_m      = E_m;
parameters.E_v      = E_v;
parameters.mu_m     = mu_m;
parameters.mu_v     = mu_v;
parameters.nu       = nu;
parameters.omega    = omega;

parameters.x_thickness = x_norm;
parameters.h_thickness = h_profile_mean;

% Model names
model_names = {'Elastic', 'Viscous', 'Maxwell'};

% Unique colors for Top 1 through Top 10
base_color = [0 0.4470 0.7410];

for j=1:3
    figure;
    ax = gca;
    hold(ax,'on');
    
    for i = 1:10
    
        parameters.w0 = calculate_mean_amp('Thwaites', i);
    
        switch j
            case 1
                w = elastic_beam_profile(L, x, parameters);
            case 2
                w = viscous_beam_profile(L, x, parameters);
            case 3
                w = maxwell_beam_profile(L, x, parameters);
        end
    
        % Top 1 = lightest, Top 10 = darkest
        alpha = 0.15 + 0.85*(i-1)/9;
    
        % Blend base color with white
        faded_color = alpha*base_color + (1-alpha)*[1 1 1];
    
        plot(ax, x_km, w, ...
            'LineWidth', 2, ...
            'Color', faded_color, ...
            'DisplayName', sprintf('Top %d', i));
    
        % Settling distance
        settling_distance = ...
            calculate_settling_distance( ...
                w, x, parameters.w0, sigma_z);
    
        % Only Top 1 and Top 10
        if i == 1 || i == 10  
            xline(ax, settling_distance, ':', ...
                sprintf('Top %d = %.1f km', i, settling_distance), ...
                'LineWidth', 1.5, ...
                'LabelVerticalAlignment', 'bottom', ...
                'LabelHorizontalAlignment', 'left', ...
                'HandleVisibility', 'off');
        end
    
    end
    
    % Mean GL distance
    xline(ax, L_mean/1e3, '--', ...
        sprintf('Mean GL Distance = %.1f km', L_mean/1e3), ...
        'Color', 'red', ...
        'LabelVerticalAlignment', 'middle', ...
        'LabelHorizontalAlignment', 'right', ...
        'HandleVisibility', 'off');
    
    xlabel(ax, 'Distance from grounding line (km)');
    ylabel(ax, 'Vertical displacement (m)');
    
    title(ax, sprintf('%s Beam Response', model_names{j}));
    
    legend(ax, 'Location', 'southeast');
    
    grid(ax,'on');
    box(ax,'on');
    
    filename = sprintf('%s_varying_amps.jpg', lower(model_names{j}));
    exportgraphics(gcf, fullfile(figure_dir, filename), ...
    'Resolution', 300);
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

%% Local functions for maxwell beam
function u_maxwell_peak = maxwell_beam_profile(L, x, parameters)
%%outputs u_maxwell profile for varying bed thickness
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
        @(x) maxwell_initial_guess(x, L, parameters.w0));
    
    % Solve
    options = bvpset( ...
        'RelTol', 1e-7, ...
        'AbsTol', 1e-9, ...
        'NMax', 20000, ...
        'Stats', 'on');
    
    sol = bvp4c( ...
        @(x,y) maxwell_beam_ode(x, y, parameters), ...
        @maxwell_beam_bc, ...
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
    
    u_maxwell_peak = ...
        us .* sin(parameters.omega*t_peak) + ...
        uc .* cos(parameters.omega*t_peak);
end 

function dydx = maxwell_beam_ode(x, y, p)
% Maxwell model for beam.
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
    omega = p.omega;
    % coefficients
    C1 = 1/(4*I_local*omega*p.mu_m);
    C2 = 3/(4*I_local*p.E_m);

    % Return derivatives of all eight state variables
    dydx = [
        alpha_s                         % dus/dx
        alpha_c                         % duc/dx
        C1*Mc + C2*Ms                   % dalpha_s/dx
        C2*Mc - C1*Ms                   % dalpha_c/dx
        p.rho_w*p.g*(p.w0 - us)         % dQs/dx
       -p.rho_w*p.g*uc                  % dQc/dx
        Qs                              % dMs/dx
        Qc                              % dMc/dx
    ];
end

function residual = maxwell_beam_bc(ya, yb)

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

function y_guess = maxwell_initial_guess(x, L, a)

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