%% Steady-state tidal flexure of a viscous beam
%
% Governing equation:
%
% Dv * d/dt(d^4w/dx^4) + rho_w*g*w
%     = rho_w*g*w0*sin(omega*t)
%
% Steady-periodic solution:
%
% w(x,t) = imag(W(x)*exp(i*omega*t))
%
% where W(x) is the complex spatial response.

clear;
clc;
close all;
%%
% extract key params
[h_shelf_mean,h_shelf_eff,h_front_mean,h_front_eff, h_calving_front] = extract_shelf_thickness('Thwaites Glacier');
mean_amp = calculate_mean_amp('Thwaites', 'top2');

% parameters 
rho = 1025; % kg/m^3
g = 9.81; % m/s^2
h = h_calving_front; % use the mean effective thickness of calving front
w0 = mean_amp; % m
sigma_z = 0.0054;  % m (vertical displacement uncertainty)

% parameters for elastic beam bending
I = h^3/12;
mu = 0.325; % poisson's ratio
E = 9.33e9; % Pa
D = (E*I)/(1-mu^2); % flexural rigidity
lambda = (rho*g*(1-mu^2)/(4*E*I))^(1/4);
lf = 1/lambda; % calculated flexural length

% parameters for viscous beam bending
T_hours = 12;            % Tidal period (hours), M2 example
T = T_hours*3600;           % Tidal period (s)
omega = 2*pi/T;             % Angular frequency (rad/s)

A_3 = 3.5e-25; % Pa^-3 s^-1
tau_e_ref = 0.1e6;  
eta_newt = 1/(2*A_3*tau_e_ref^2);
D_v = eta_newt*h^3/3;

% parameters for non-Newtonian FEM
R_gas_local = 8.314; % Gas constant
eps_bg = 1e-10; % background extensional strain rate, 1/s
kappa_bg = eps_bg/(h/2); % equivalent background curvature rate
T_ice = 263.0;  

% flow law parameters
A_1 = 1/(2*eta_newt); % newtonian flow parameter
%% Profile specifications
n_lf = 30;

Nel = 160; % number of elements
L = n_lf*lf; % length of ice shelf 
x = linspace(0,L,Nel);
x_km = x/1000;
%% solving analytical viscous profile

% calculate beta 
beta = (rho*g/(omega*D_v))^(1/4);

% Characteristic viscous flexural length
ell_v = 1/beta;

% characteristic roots
lambda_all = beta*exp(1i*(pi/8 + (0:3)*pi/2));
lambda_decay = lambda_all(real(lambda_all) < 0);
lambda1 = lambda_decay(1);
lambda2 = lambda_decay(2);

% calculates constants
C1 =  w0*lambda2/(lambda1-lambda2);
C2 = -w0*lambda1/(lambda1-lambda2);

% complex spatial response
W = w0 ...
    + C1.*exp(lambda1.*x) ...
    + C2.*exp(lambda2.*x);

% extract local amplitude and phase
amplitude = abs(W);                 % Local tidal amplitude (m)

% W = amplitude*exp(-i*phase_lag)
phase_lag = -angle(W);              % Phase lag (rad)
phase_lag = unwrap(phase_lag);      % Remove artificial 2*pi jumps
phase_lag_deg = rad2deg(phase_lag); % Phase lag (degrees)

% The phase is not meaningful exactly where amplitude is zero.
phase_lag_deg(amplitude < 1e-8*w0) = NaN;

% Peak high tide occurs at t = T/4
t_peak = T/4;

% Instantaneous physical displacement profile
w_analytical = imag(W .* exp(1i*omega*t_peak));

%% Power Law viscous tidal flexure model
% ================================================================
%  Flow-law cases
% ================================================================
n_values = [1.0];
A_values = [A_1];

flow_labels = { ...
    'n = 1 (Newtonian)'};

for j = 1:numel(n_values)

    n_val = n_values(j);
    A_val = A_values(j);

    eps_dot_ref = A_val*tau_e_ref^n_val;

    eta_eff_ref = ...
        1/(2*A_val*tau_e_ref^(n_val - 1));

    fprintf(['  %s: A = %.2e, ' ...
             'eps_dot(0.1 MPa) = %.2e 1/s, ' ...
             'eta = %.2e Pa s\n'], ...
        flow_labels{j}, A_val, eps_dot_ref, eta_eff_ref);

end

% ================================================================
% Spatial grid
% ================================================================

he = L/Nel;        % Element length, m

x_nodes = linspace(0, L, Nel + 1)';

% Each node has:
%   1. vertical velocity
%   2. slope rate
ndof = 2*(Nel + 1);


% ================================================================
% Two-point Gauss quadrature
% ================================================================

gp = [ ...
    0.5 - sqrt(3)/6, ...
    0.5 + sqrt(3)/6];

gw = [0.5, 0.5];


% Hermite shape functions

N_gp = zeros(2,4);
d2N_gp = zeros(2,4);

for q = 1:2

    xi = gp(q);

    % Cubic Hermite shape functions
    N_gp(q,:) = [ ...
        1 - 3*xi^2 + 2*xi^3, ...
        xi*(1 - xi)^2, ...
        3*xi^2 - 2*xi^3, ...
        xi^2*(xi - 1)];

    % Second derivatives with respect to xi
    d2N_gp(q,:) = [ ...
        -6 + 12*xi, ...
        -4 + 6*xi, ...
         6 - 12*xi, ...
        -2 + 6*xi];

end


% Convert shape functions to physical coordinates

% The slope degrees of freedom require factors of element length.
scale = [1, he, 1, he];

N_phys = N_gp .* scale;

% d^2/dx^2 = (1/he^2) d^2/dxi^2
B_phys = (d2N_gp/he^2) .* scale;


% ================================================================
% Reference element matrices
% ================================================================

% K_ref(:,:,q) contains:
%
%     weight * he * B' * B
%
% for Gauss point q.

K_ref = zeros(4,4,2);

% F_ref(:,q) contains:
%
%     weight * he * N'
%
F_ref = zeros(4,2);

% Hydrostatic restoring matrix:
%
%     integral N' N dx

M_ref = zeros(4,4);

for q = 1:2

    Bq = B_phys(q,:);
    Nq = N_phys(q,:);

    K_ref(:,:,q) = ...
        gw(q)*he*(Bq'*Bq);

    F_ref(:,q) = ...
        gw(q)*he*Nq';

    M_ref = M_ref + ...
        gw(q)*he*(Nq'*Nq);

end


% ================================================================
% Element-to-global degree-of-freedom mapping
% ================================================================

element_dofs = zeros(Nel,4);

% These arrays are used for efficient sparse-matrix assembly.
row_indices = zeros(Nel,4,4);
col_indices = zeros(Nel,4,4);

for e = 1:Nel

    % Element e connects node e and node e+1.
    %
    % MATLAB degree-of-freedom ordering:
    %
    % node 1: [w_dot_1, slope_dot_1]
    % node 2: [w_dot_2, slope_dot_2]
    %
    dofs = [ ...
        2*e - 1, ...
        2*e, ...
        2*e + 1, ...
        2*e + 2];

    element_dofs(e,:) = dofs;

    row_indices(e,:,:) = repmat(dofs',1,4);
    col_indices(e,:,:) = repmat(dofs,4,1);

end

row_indices_flat = row_indices(:);
col_indices_flat = col_indices(:);

% Store all model information in one structure
model.h = h;
model.w0 = w0;
model.T = T;
model.omega = omega;

model.rho = rho;
model.g = g;

model.Nel = Nel;
model.he = he;
model.ndof = ndof;
model.x_nodes = x_nodes;

model.gp = gp;
model.gw = gw;

model.N_phys = N_phys;
model.B_phys = B_phys;

model.K_ref = K_ref;
model.F_ref = F_ref;
model.M_ref = M_ref;

model.element_dofs = element_dofs;
model.row_indices_flat = row_indices_flat;
model.col_indices_flat = col_indices_flat;

model.kappa_bg = kappa_bg;
model.D_v = D_v;


% ================================================================
% Version 1: Physical flow-law values
% ================================================================

fprintf('\n--- Version 1: Physical A values ---\n');

nl_phys = struct();

for j = 1:numel(n_values)

    n_val = n_values(j);
    A_val = A_values(j);

    D_n = compute_D_n(A_val, n_val, h);

    result = run_power_law( ...
        n_val, ...
        D_n, ...
        flow_labels{j}, ...
        6, ...
        model);

    field_name = matlab.lang.makeValidName( ...
        sprintf('n_%g',n_val));

    nl_phys.(field_name) = result;

end


% ================================================================
% Version 2: Match D_eff(kappa_bg) to D_v
% ================================================================

fprintf('\n--- Version 2: Matched D_eff(kappa_bg) = D_v ---\n');

nl_matched = struct();

for j = 1:numel(n_values)

    n_val = n_values(j);

    % D_eff = D_n*kappa_bg^(1/n - 1)
    %
    % Setting D_eff(kappa_bg) = D_v gives:
    %
    % D_n = D_v*kappa_bg^(1 - 1/n)

    D_n_matched = ...
        D_v*kappa_bg^(1 - 1/n_val);

    matched_label = ...
        sprintf('%s, matched',flow_labels{j});

    result = run_power_law( ...
        n_val, ...
        D_n_matched, ...
        matched_label, ...
        6, ...
        model);

    field_name = matlab.lang.makeValidName( ...
        sprintf('n_%g',n_val));

    nl_matched.(field_name) = result;

end

fprintf('\nPower-law viscous simulations complete.\n');

% Extract results
x_powerlaw_km = x_nodes/1000;

w_n1 = nl_phys.n_1.high_tide_profile;

% Make sure all profiles are column vectors
w_n1 = w_n1(:); 

%% Shooting Method
% Parameters
parameters.rho_w = rho;
parameters.g     = g;
parameters.a     = w0;
parameters.I     = I;
parameters.mu_m  = eta_newt;
parameters.omega = omega;

% Initial guess
x_mesh = linspace(0, L, 300);

solinit = bvpinit( ...
    x_mesh, ...
    @(x_temp) viscous_initial_guess(x_temp, L, w0));

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

x_shoot = linspace(0, L, 1000);
Y = deval(sol, x_shoot);
x_shoot_km = x_shoot/1e3;

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

t_peak = T/4;

u_shooting_peak = ...
    us .* sin(omega*t_peak) + ...
    uc .* cos(omega*t_peak);
%% Settling-distance calculations

% calculate settling distance
viscous_settling_distance = calculate_settling_distance(w_analytical, x, w0, sigma_z);

% Power-law FEM profiles
settling_distance_n1 = ...
    calculate_settling_distance( ...
        w_n1, x_nodes, w0, sigma_z);

settling_distance_shooting = calculate_settling_distance(u_shooting_peak, x_shoot, w0, sigma_z);
%% Plotting Profile comparisons betweebn analytical, FEM, and Shooting method

figure;
hold on;

% viscous
plot(x_km, w_analytical, ...
    'LineStyle', '-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Analytical'); 

% Power-law FEM profiles
plot(x_powerlaw_km, w_n1, ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'DisplayName', 'FEM');

% shooting method profile
plot(x_shoot_km, u_shooting_peak, ...
    'LineStyle', ':', ...
    'LineWidth', 2, ...
    'DisplayName', 'Shooting'); 

% Far-field steady-state tidal displacement
yline(w0, '--', '$w_0$', ...
    'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left');

% Upper and lower uncertainty bounds around steady state
upper_bound = w0 + sigma_z;
lower_bound = w0 - sigma_z;

yline(upper_bound, ':', ...
    '$\pm\sigma_z$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right');
yline(lower_bound, ':');

% Settling-distance markers
xline(viscous_settling_distance, ':', ...
    sprintf('Analytical = %.1f km', viscous_settling_distance), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left');

xline(settling_distance_n1, ':', ...
    sprintf('FEM = %.1f km', settling_distance_n1), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left');

xline(settling_distance_shooting, ':', ...
    sprintf('Shooting = %.1f km', settling_distance_shooting), ...
    'LabelVerticalAlignment', 'top', ...
    'LabelHorizontalAlignment', 'left');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Deflection at Peak Tide Thwaites Glacier');

xlim([0, L/1000]);

% Axes and labels

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Viscous Peak Displacement Profile: Analytical, FEM, Shooting Method');

xlim([0, max([x_km(:); x_powerlaw_km(:)])]);

% Set y-limits using every plotted profile

all_y_values = [ ...
    w_analytical(:); ...
    w_n1(:);...
    upper_bound; ...
    lower_bound];

w_range = max(all_y_values) - min(all_y_values);

if w_range == 0
    vertical_padding = ...
        max(abs(all_y_values))*0.1;
else
    vertical_padding = ...
        0.20*w_range;
end

ylim([ ...
    min(all_y_values) - vertical_padding, ...
    max(all_y_values) + vertical_padding]);


% Plot formatting

legend('Analytical', 'FEM', 'Shooting',...
    'Location', 'southeast');

box on;

set(gca, ...
    'FontSize', 12, ...
    'LineWidth', 1);

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/viscous_comparisons';
exportgraphics(gcf, fullfile(figure_dir,'viscous_method_comparison.jpg'), ...
'Resolution',300);

%% Plotting Peak tide displacement profile

% Peak high tide occurs at t = T/4
t_peak = T/4;

% Instantaneous physical displacement profile
w_peak = imag(W .* exp(1i*omega*t_peak));

% Settling-distance calculations
viscous_settling_distance = calculate_settling_distance(w_peak, x, w0, sigma_z);

figure;
hold on;
plot(x_km, w_peak, 'LineWidth', 2); % viscous

% Far-field steady-state tidal displacement
yline(w0, '--', '$w_0$', ...
    'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left');

% Upper and lower uncertainty bounds around steady state
upper_bound = w0 + sigma_z;
lower_bound = w0 - sigma_z;

yline(upper_bound, ':', ...
    '$\pm\sigma_z$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right');
yline(lower_bound, ':');

% Settling-distance markers
xline(viscous_settling_distance, ':', ...
    sprintf('Viscous Settling distance = %.1f km', viscous_settling_distance), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Deflection at Peak Tide Thwaites Glacier');

xlim([0, L/1000]);

% Include the displacement profile and uncertainty bounds in the y-limits
all_y_values = [w_peak(:); upper_bound; lower_bound];

w_range = max(all_y_values) - min(all_y_values);

% Avoid zero padding if the plotted values are nearly constant
if w_range == 0
    vertical_padding = max(abs(all_y_values))*0.1;
else
    vertical_padding = 0.20*w_range;
end

ylim([min(all_y_values) - vertical_padding, ...
      max(all_y_values) + vertical_padding]);

grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);

% figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending';
% exportgraphics(gcf, fullfile(figure_dir,'elastic_viscous_beam_bending_thwaites.jpg'), ...
% 'Resolution',300);

%% Plot time series for one cycle
figure;
hold on;

% Tidal phases to plot
nPhases = 8;

% Plot 0°, 45°, ..., 315°. Do not plot 360° because it duplicates 0°.
phases_deg = linspace(0, 360, nPhases + 1);
phases_deg(end) = [];


% Generate a 256-color twilight colormap
cmap = slanCM('twilight_s',256);
nColors = size(cmap,1);

% Apply to current axes
colormap(gca, cmap);

% Colorbar spans 0–360°
clim([0 360]);


% Store all profiles for determining y-limits
all_y_values = [];

for j = 1:nPhases

    phase_deg = phases_deg(j);

    % Convert phase angle to time
    t = (phase_deg/360)*T;

    % Instantaneous physical displacement profile
    w = imag(W .* exp(1i*omega*t));

    % Select the line color from the continuous colormap
    color_position = phase_deg/360;
    color_index = 1 + round(color_position*(nColors - 1));
    line_color = cmap(color_index,:);

    plot(x_km, w, ...
        'Color', line_color, ...
        'LineWidth', 2);

    all_y_values = [all_y_values; w(:)];

end

% Axis labels
xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Thwaites Glacier Tidal Deflection');

xlim([0, L/1000]);

% Set y-limits from all phases
w_min = min(all_y_values);
w_max = max(all_y_values);
w_range = w_max - w_min;

if w_range == 0
    vertical_padding = 0.1*max(abs(all_y_values));

    if vertical_padding == 0
        vertical_padding = 1;
    end
else
    vertical_padding = 0.20*w_range;
end

ylim([w_min - vertical_padding, ...
      w_max + vertical_padding]);

% Continuous phase colorbar
cb = colorbar;
cb.Label.String = 'Tidal phase (degrees)';
cb.Ticks = 0:90:360;
cb.TickLabels = {'0^\circ','90^\circ','180^\circ', ...
                 '270^\circ','360^\circ'};

grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);

% 
% figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending';
% exportgraphics(gcf, fullfile(figure_dir,'viscous_tidal_phases.jpg'), ...
% 'Resolution',300);

%% Creating animation to show phases
% Number of animation frames
nFrames = 120;

% Go through one complete tidal cycle
phases_deg = linspace(0, 360, nFrames);

% Twilight colormap
cmap = slanCM('twilight_s',256);
nColors = size(cmap,1);

% -------------------------------------------------
% Determine fixed y limits before animation
% -------------------------------------------------

% Maximum possible displacement envelope
w_amp = abs(W);

w_max = max(w_amp);

vertical_padding = 0.15*w_max;

y_limits = [ ...
    -w_max - vertical_padding, ...
     w_max + vertical_padding];

% -------------------------------------------------
% Create figure
% -------------------------------------------------

fig = figure( ...
    'Position',[100 100 1000 650], ...
    'Color','w');

ax = axes(fig);
hold(ax,'on');

% Initial phase
phase_deg = phases_deg(1);
t_current = (phase_deg/360)*T;

w = imag(W .* exp(1i*omega*t_current));

% Initial color
color_index = 1;
line_color = cmap(color_index,:);

% Plot displacement
h = plot(x_km, w, ...
    'Color',line_color, ...
    'LineWidth',3);

% Grounding line / zero displacement
yline(0,'k--', ...
    'HandleVisibility','off');

% Formatting
xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');

title('Thwaites Glacier Tidal Deflection');

xlim([0 L/1000]);
ylim(y_limits);

grid on;
box on;

set(gca, ...
    'FontSize',12, ...
    'LineWidth',1);

% Colorbar
colormap(gca,cmap);
clim([0 360]);

cb = colorbar;
cb.Label.String = 'Tidal phase (degrees)';

cb.Ticks = 0:90:360;
cb.TickLabels = { ...
    '0^\circ', ...
    '90^\circ', ...
    '180^\circ', ...
    '270^\circ', ...
    '360^\circ'};

% Moving phase label
phase_text = text( ...
    0.03,0.92, ...
    'Phase = 0^\circ', ...
    'Units','normalized', ...
    'FontSize',14, ...
    'FontWeight','bold');

% -------------------------------------------------
% GIF settings
%%-------------------------------------------------

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/animations';

if ~exist(figure_dir,'dir')
    mkdir(figure_dir);
end

fname = fullfile(figure_dir,'viscous_phases');

framerate = 25;

% -------------------------------------------------
% Animate
% -------------------------------------------------

for j = 1:nFrames

    % Current tidal phase
    phase_deg = phases_deg(j);

    % Convert phase to time
    t_current = (phase_deg/360)*T;

    % Calculate instantaneous displacement
    w = imag(W .* exp(1i*omega*t_current));

    % Select color corresponding to phase
    color_position = phase_deg/360;

    color_index = 1 + ...
        round(color_position*(nColors-1));

    line_color = cmap(color_index,:);

    % Update displacement
    h.YData = w;

    % Update line color
    h.Color = line_color;

    % Update phase label
    phase_text.String = ...
        sprintf('Phase = %.0f^\\circ',phase_deg);

    drawnow limitrate

    % Write GIF
    if j == 1
        gifanim(fname,framerate,fig,1);
    else
        gifanim(fname,framerate,fig,0);
    end

end

%% varying the tau_e valeus

tau_e_array = [0.05e6, 0.1e6, 0.15e6, 0.2e6, 0.25e6, 0.3e6]; % effective stress
figure;
hold on; 

for i=1:length(tau_e_array)
    tau_e = tau_e_array(i);
    eta = 1 / (2 * A_3 * tau_e^2); % Update viscosity for current tau_e
    Dv = eta * h^3 / 3;          % Calculate viscous flexural rigidity

    % calculate beta 
    beta = (rho*g/(omega*Dv))^(1/4);
    
    % Characteristic viscous flexural length
    ell_v = 1/beta;
    
    % characteristic roots
    lambda_all = beta*exp(1i*(pi/8 + (0:3)*pi/2));
    lambda_decay = lambda_all(real(lambda_all) < 0);
    lambda1 = lambda_decay(1);
    lambda2 = lambda_decay(2);
    
    % calculates constants
    C1 =  w0*lambda2/(lambda1-lambda2);
    C2 = -w0*lambda1/(lambda1-lambda2);
    
    % complex spatial response
    W = w0 ...
        + C1.*exp(lambda1.*x) ...
        + C2.*exp(lambda2.*x);
    
    % extract local amplitude and phase
    amplitude = abs(W);                 % Local tidal amplitude (m)
    
    % W = amplitude*exp(-i*phase_lag)
    phase_lag = -angle(W);              % Phase lag (rad)
    phase_lag = unwrap(phase_lag);      % Remove artificial 2*pi jumps
    phase_lag_deg = rad2deg(phase_lag); % Phase lag (degrees)
    
    % The phase is not meaningful exactly where amplitude is zero.
    phase_lag_deg(amplitude < 1e-8*w0) = NaN;

    % Peak high tide occurs at t = T/4
    t_peak = T/4;
    
    % Instantaneous physical displacement profile
    w_peak = imag(W .* exp(1i*omega*t_peak));
    
    % Settling-distance calculations
    viscous_settling_distance = calculate_settling_distance(w_peak, x, w0, sigma_z);

    plot(x_km, w_peak, ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('\\tau_e = %.2fMPa\', tau_e/1e6));
    % Settling-distance markers
    xline(viscous_settling_distance, ':', ...
        sprintf('\\tau_e = %.2fMPa: %.1fkm\', tau_e/1e6, viscous_settling_distance),...
        'LabelVerticalAlignment', 'middle', ...
        'LabelHorizontalAlignment', 'left', ...
        'HandleVisibility', 'off');
end


% Far-field steady-state tidal displacement
yline(w0, '--', '$w_0$', ...
    'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

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

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Viscous Deflection at Peak Tide Thwaites Glacier');

xlim([0, L/1000]);

% Include the displacement profile and uncertainty bounds in the y-limits
all_y_values = [w_peak(:); upper_bound; lower_bound];

w_range = max(all_y_values) - min(all_y_values);

% Avoid zero padding if the plotted values are nearly constant
if w_range == 0
    vertical_padding = max(abs(all_y_values))*0.1;
else
    vertical_padding = 0.20*w_range;
end

ylim([min(all_y_values) - vertical_padding, ...
      max(all_y_values) + vertical_padding]);

legend('Location','southeast')

grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/param_variations';
exportgraphics(gcf, fullfile(figure_dir,'varied_tau_viscous.jpg'), ...
'Resolution',300);


%% Functions for FEM solver

function A_Glen = A_Glen_from_GK( ...
    A0_MPa, n, Q, T, d, p, R_gas)
% A_GLEN_FROM_GK
% Apply temperature, grain-size, and stress-unit conversions.
%
% A0_MPa:
%     Pre-exponential factor with stress measured in MPa.
%
% d:
%     Grain size in meters. Pass [] when grain size is not used.
%
% Important:
% The final GK-to-Glen multiplicative factor should be independently
% verified against the tensor conventions used in the model.

    % Arrhenius temperature dependence
    A_GK_MPa = ...
        A0_MPa*exp(-Q/(R_gas*T));

    % Grain-size dependence
    if ~isempty(d) && p > 0
        A_GK_MPa = A_GK_MPa/d^p;
    end

    % Convert MPa^-n to Pa^-n
    A_GK_Pa = ...
        A_GK_MPa/(1e6)^n;

    % Convention conversion used in the Python code
    A_Glen = ...
        3^((n + 1)/2)/2*A_GK_Pa;

end

function D_n = compute_D_n(A_Glen, n_nl, h_plate)
% COMPUTE_D_N
% Nonlinear viscous bending coefficient for the relation
%
% M = D_n*|kappa_dot|^(1/n - 1)*kappa_dot

    D_n = ...
        (4*n_nl/(1 + 2*n_nl)) ...
        *(h_plate/2)^(1/n_nl + 2) ...
        /A_Glen^(1/n_nl);

end



function wdot_vec = solve_wdot_fem( ...
    w_curr, ...
    w_tide_scalar, ...
    D_n, ...
    n_nl, ...
    dt, ...
    wdot_previous, ...
    picard_iterations, ...
    model)
% SOLVE_WDOT_FEM
% Solve for nodal vertical velocities and nodal slope rates during one
% time step.
%
% Global unknown ordering:
%
% [w_dot_1;
%  slope_dot_1;
%  w_dot_2;
%  slope_dot_2;
%  ... ]

    Nel = model.Nel;
    ndof = model.ndof;

    B_phys = model.B_phys;
    K_ref = model.K_ref;
    F_ref = model.F_ref;
    M_ref = model.M_ref;

    rho = model.rho;
    g = model.g;

    gp = model.gp;

    kappa_bg = model.kappa_bg;

    element_dofs = model.element_dofs;

    row_indices_flat = model.row_indices_flat;
    col_indices_flat = model.col_indices_flat;


    % Initial Picard guess

    if isempty(wdot_previous)
        wdot_vec = zeros(ndof,1);
    else
        wdot_vec = wdot_previous;
    end


    % Linear case requires only one matrix solve

    if n_nl > 1
        number_iterations = picard_iterations;
    else
        number_iterations = 1;
    end


    % Picard nonlinear iteration

    for iteration = 1:number_iterations

        % One effective bending coefficient at each Gauss point
        % of every element.
        D_eff = D_n*ones(Nel,2);


        % Calculate nonlinear effective bending resistance

        if n_nl > 1

            for e = 1:Nel

                dofs = element_dofs(e,:);

                % Element velocity degrees of freedom:
                %
                % [left velocity;
                %  left slope rate;
                %  right velocity;
                %  right slope rate]
                dv = wdot_vec(dofs);

                for q = 1:2

                    % Curvature rate:
                    %
                    % kappa_dot = d^2(w_dot)/dx^2
                    kappa_dot = ...
                        B_phys(q,:)*dv;

                    % Smooth regularization
                    kappa_dot_eff = ...
                        sqrt(kappa_dot^2 + kappa_bg^2);

                    % Nonlinear effective bending coefficient
                    D_eff(e,q) = ...
                        D_n ...
                        *kappa_dot_eff^(1/n_nl - 1);

                end

            end

        end


        % Assemble global matrix element by element

        K_global = sparse(ndof,ndof);
        
        for e = 1:Nel
        
            % Four global degrees of freedom for this element
            dofs = element_dofs(e,:);
        
            % Local 4-by-4 element matrix
            K_e = ...
                  D_eff(e,1)*K_ref(:,:,1) ...
                + D_eff(e,2)*K_ref(:,:,2) ...
                + dt*rho*g*M_ref;
        
            % Add local matrix into the global matrix
            K_global(dofs,dofs) = ...
                K_global(dofs,dofs) + K_e;
        
        end


        % Current displacement at Gauss points

        w_gp = zeros(Nel,2);

        for q = 1:2

            xi = gp(q);

            % This reproduces the Python code exactly:
            % linear interpolation between nodal displacements.
            w_gp(:,q) = ...
                (1 - xi)*w_curr(1:end-1) ...
                + xi*w_curr(2:end);

        end


        % Hydrostatic load

        load_gp = ...
            rho*g*(w_tide_scalar - w_gp);


        % Element force vectors

        F_element = zeros(Nel,4);

        for e = 1:Nel

            F_element(e,:) = ...
                (load_gp(e,1)*F_ref(:,1) ...
                + load_gp(e,2)*F_ref(:,2))';

        end


        % Assemble global force vector

        F_global = accumarray( ...
            element_dofs(:), ...
            F_element(:), ...
            [ndof,1], ...
            @sum, ...
            0);


        % Clamped grounding-line boundary conditions

        constrained_dofs = [1,2];

        % Degree of freedom 1:
        %     w_dot(0) = 0
        %
        % Degree of freedom 2:
        %     d(w_dot)/dx at x=0 = 0

        for bc = constrained_dofs

            K_global(bc,:) = 0;
            K_global(:,bc) = 0;

            K_global(bc,bc) = 1;
            F_global(bc) = 0;

        end


        % Solve the linear system

        wdot_vec = K_global\F_global;

    end

end



function result = run_power_law( ...
    n_nl, ...
    D_n, ...
    label, ...
    number_cycles, ...
    model)
% RUN_POWER_LAW
% Run a power-law viscous simulation for multiple tidal cycles.

    kappa_bg = model.kappa_bg;
    D_v = model.D_v;

    Nel = model.Nel;
    ndof = model.ndof;

    T = model.T;
    omega = model.omega;
    w0 = model.w0;

    he = model.he;


    % Background effective rigidity

    D_eff_background = ...
        D_n*kappa_bg^(1/n_nl - 1);

    fprintf(['  %s: D_n = %.2e, ' ...
             'D_eff(kappa_bg) = %.2e, ' ...
             'ratio to D_v = %.2f\n'], ...
        label, ...
        D_n, ...
        D_eff_background, ...
        D_eff_background/D_v);


    % Time discretization

    dt = 30.0;

    steps_per_cycle = ...
        round(T/dt);

    number_steps = ...
        number_cycles*steps_per_cycle;

    % Save approximately 200 profiles during the last cycle.
    save_every = ...
        max(1,floor(steps_per_cycle/200));

    save_start = ...
        (number_cycles - 1)*steps_per_cycle + 1;


    % Initial conditions

    % Only nodal displacement is explicitly stored here.
    w = zeros(Nel + 1,1);

    wdot_previous = [];

    maximum_saved_profiles = ...
        ceil(steps_per_cycle/save_every) + 1;

    w_history = ...
        zeros(maximum_saved_profiles,Nel + 1);

    time_history = ...
        zeros(maximum_saved_profiles,1);

    save_count = 0;


    % Time stepping

    for step = 1:number_steps

        % MATLAB indexing begins at 1, but physical time begins at zero.
        current_time = ...
            (step - 1)*dt;

        w_tide = ...
            w0*sin(omega*current_time);


        % Solve for vertical velocity

        wdot_vec = solve_wdot_fem( ...
            w, ...
            w_tide, ...
            D_n, ...
            n_nl, ...
            dt, ...
            wdot_previous, ...
            3, ...
            model);


        % Save velocity as next initial Picard guess

        wdot_previous = wdot_vec;


        % Update nodal displacement

        % Odd MATLAB indices 1,3,5,... are vertical velocities.
        vertical_velocity = ...
            wdot_vec(1:2:ndof);

        w = ...
            w + dt*vertical_velocity;

        % Explicitly enforce grounding-line displacement.
        w(1) = 0;


        % Save only the final cycle

        if step >= save_start && ...
                mod(step - save_start,save_every) == 0

            save_count = save_count + 1;

            w_history(save_count,:) = w';
            time_history(save_count) = current_time;

        end

    end


    % Remove unused preallocated rows

    w_history = ...
        w_history(1:save_count,:);

    time_history = ...
        time_history(1:save_count);


    % Find profile nearest positive high tide

    phase = ...
        omega*time_history;

    [~,high_tide_index] = ...
        max(sin(phase));

    high_tide_profile = ...
        w_history(high_tide_index,:);


    % Diagnostics

    maximum_displacement = ...
        max(abs(w_history),[],'all');

    offshore_displacement = ...
        high_tide_profile(end);

    offshore_slope = ...
        (high_tide_profile(end) ...
        - high_tide_profile(end - 1))/he;

    fprintf(['    max|w| = %.3f m, ' ...
             'w(L) = %.4f m, ' ...
             'dw/dx(L) = %.2e\n'], ...
        maximum_displacement, ...
        offshore_displacement, ...
        offshore_slope);

    % Return results

    result.n = n_nl;
    result.D_n = D_n;

    result.w_history = w_history;
    result.time_history = time_history;

    result.final_displacement = w;
    result.final_wdot = wdot_previous;

    result.high_tide_profile = ...
        high_tide_profile';

end

%% Local functions for viscous beam 
function dydx = viscous_beam_ode(~, y, p)
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

    % Material coefficients
    omega = p.omega;
    I = p.I; 
    mu_m = p.mu_m; 
    Dv = 4*mu_m*I;


    % Return derivatives of all eight state variables
    dydx = [
        alpha_s                         % dus/dx
        alpha_c                         % duc/dx
        Mc/(Dv*omega)                     % dalpha_s/dx
        -Ms/(Dv*omega)                     % dalpha_c/dx
        p.rho_w*p.g*(p.a - us)          % dQs/dx
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