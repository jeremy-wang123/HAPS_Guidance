%% Redoing Reeh et al., 2003 Analysis for viscoelastic flexure
clear;clc;
figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/reeh';

%% Parameter values
[h_shelf_mean,h_shelf_eff,h_front_mean,h_front_eff, h_calving_front] = extract_shelf_thickness('Thwaites Glacier');
mean_amp = calculate_mean_amp('Thwaites', 'top2');

rho_w = 1025; % kg/m^3
g = 9.81; % m/s^2
h = h_calving_front; % use the mean effective thickness of calving front
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

% flexural length
I_analytical = h^3 / 12;          % second moment per unit width, m^3
D_analytical = E_m*I_analytical / (1 - nu^2);  % flexural rigidity per unit width, N m
lambda = (rho_w*g/(4*D_analytical))^(1/4); 
L = 30/lambda;

%% Comparing numerical method and analytical solution to elastic beam
x = linspace(0, L, 1000);

% Flexural properties
I = h^3 / 12;          % second moment per unit width, m^3
D = E_m*I / (1 - nu^2);  % flexural rigidity per unit width, N m

% Characteristic elastic decay parameter
lambda = (rho_w*g/(4*D))^(1/4);
w_analytical = w0 .* ...
    (1 - exp(-lambda*x) .* ...
    (cos(lambda*x) + sin(lambda*x)));

% Initial mesh and initial guess
x_mesh = linspace(0, L, 300);

solinit = bvpinit( ...
    x_mesh, ...
    @(x) initial_guess(x, L, w0));

% Solve boundary-value problem
parameters.rho_w = rho_w;
parameters.g     = g;
parameters.D     = D;
parameters.w0    = w0;

options = bvpset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-9, ...
    'Stats', 'on');

sol = bvp4c( ...
    @(x,y) beam_ode(x, y, parameters), ...
    @beam_boundary_conditions, ...
    solinit, ...
    options);

% Evaluate numerical solution
x = linspace(0, L, 1000);

Y = deval(sol, x);

w_elastic     = Y(1,:);
theta = Y(2,:);
M     = Y(3,:);
Q     = Y(4,:);

x_km = x/1e3;
% Plot displacement for elastic
figure;
hold on;
grid on;
box on;

plot(x_km, w_elastic, ...
    'LineWidth', 2, ...
    'DisplayName', 'Finite numerical beam');

plot(x_km, w_analytical, '--', ...
    'LineWidth', 2, ...
    'DisplayName', 'Semi-infinite analytical beam');

yline(w0, ':', ...
    '$w_0$', ...
    'Interpreter', 'latex', ...
    'HandleVisibility', 'off');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title(sprintf('Uniform thickness: h = %.0f m, L = %.0f km', ...
    h, L/1e3));

legend('Location', 'southeast');

%% Reeh Viscoelastic Beam
% Beam geometry
I = h^3/12;

% Parameters
parameters.rho_w = rho_w;
parameters.g     = g;
parameters.a     = w0;
parameters.I     = I;

parameters.E_m   = E_m;
parameters.E_v   = E_v;
parameters.mu_v  = mu_v;
parameters.mu_m  = mu_m;

parameters.omega = omega;

% Initial guess
x_mesh = linspace(0, L, 300);

solinit = bvpinit( ...
    x_mesh, ...
    @(x) viscoelastic_initial_guess(x, L, w0));

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

x = linspace(0, L, 1000);
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

t_peak = T/4;

u_viscoelastic_peak = ...
    us .* sin(omega*t_peak) + ...
    uc .* cos(omega*t_peak);

%% Viscous response
% Beam geometry
I = h^3/12;

% Parameters
parameters.rho_w = rho_w;
parameters.g     = g;
parameters.a     = w0;
parameters.I     = I;
parameters.mu_m  = mu_m;
parameters.omega = omega;

% Initial guess
x_mesh = linspace(0, L, 300);

solinit = bvpinit( ...
    x_mesh, ...
    @(x) viscous_initial_guess(x, L, w0));

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

x = linspace(0, L, 1000);
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

t_peak = T/4;

u_viscous_peak = ...
    us .* sin(omega*t_peak) + ...
    uc .* cos(omega*t_peak);
%% Displacement profile at peak high tide
figure;
hold on;
grid on;
box on;

plot(x_km, w_elastic, ...
    'LineWidth', 2, ...
    'DisplayName', 'Elastic beam');

plot(x_km, u_viscoelastic_peak, ...
    'LineWidth', 2, ...
    'DisplayName', 'Viscoelastic beam');

plot(x_km, u_viscous_peak, ...
    'LineWidth', 2, ...
    'DisplayName', 'Viscous beam');

yline(w0, '--', ...
    '$w_0$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','left',...
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

elastic_settling = ...
    calculate_settling_distance( ...
        w_elastic, x, w0, sigma_z);

viscoelastic_settling = ...
    calculate_settling_distance( ...
        u_viscoelastic_peak, x, w0, sigma_z);

viscous_settling = ...
    calculate_settling_distance( ...
        u_viscous_peak, x, w0, sigma_z);

xline(elastic_settling, ':', ...
    sprintf('Elastic = %.1f km', elastic_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(viscoelastic_settling, ':', ...
    sprintf('Viscoelastic = %.1f km', viscoelastic_settling), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(viscous_settling, ':', ...
    sprintf('Viscous = %.1f km', viscous_settling), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Displacement Profile at Peak High Tide Numerical');
legend('Location','southeast')

exportgraphics(gcf, fullfile(figure_dir,'reeh_viscoelastic_thwaites.jpg'), ...
'Resolution',300);


%% time series plot for one cycle (viscoelastic)
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
    t = (phase_deg/360)*T
    
    % instantaneous displacement profile
    w = ...
        us .* sin(omega*t) + ...
        uc .* cos(omega*t);


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
title('Thwaites Glacier Viscoelastic Tidal Deflection');

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
exportgraphics(gcf, fullfile(figure_dir,'reeh_viscoelastic_tidal_phases.jpg'), ...
'Resolution',300);

%% Local functions for elastic beam

function dydx = beam_ode(~, y, p)
% State vector:
% y(1) = vertical displacement, w
% y(2) = tilt, theta = dw/dx
% y(3) = bending moment, M
% y(4) = transverse shear force, Q

    w     = y(1);
    theta = y(2);
    M     = y(3);
    Q     = y(4);

    dydx = [
        theta
        M / p.D
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

%% Local functions for viscoelastic beam 
function dydx = viscoelastic_beam_ode(~, y, p)
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
    C1 = 3*p.mu_v/p.E_v;

    C2 = 2.25*p.mu_v / ...
        (p.E_v*p.E_m*p.I);

    C3 = 3/(4*p.I) * ...
        (1/p.E_v + ...
         1/p.E_m + ...
         p.mu_v/(p.mu_m*p.E_v));

    C4 = 2.25/(p.mu_m*p.I);

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
        p.rho_w*p.g*(p.a - us)          % dQs/dx
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
    Dv = 4*I*mu_m;


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







