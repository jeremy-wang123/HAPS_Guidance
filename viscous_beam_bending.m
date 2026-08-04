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

%% parameter extraction
% extract key params
[h_shelf_mean,h_shelf_eff,h_front_mean,h_front_eff] = extract_shelf_thickness('Thwaites Glacier');
mean_amp = calculate_mean_amp('Thwaites', 'top2');

% parameters 
rho_w = 1025; % kg/m^3
g = 9.81; % m/s^2
h = h_front_eff; % use the mean effective thickness of calving front
w0 = mean_amp; % m
sigma_z = 0.002;  % m (vertical displacement uncertainty)

% parameters for viscous beam bending
T_hours = 12;            % Tidal period (hours), M2 example
T = T_hours*3600;           % Tidal period (s)
omega = 2*pi/T;             % Angular frequency (rad/s)

A = 3.5e-25; % Pa^-3 s^-1
tau_e = 0.1e6; % Pa (effective stress)
eta = 1 / (2*A*tau_e^2); % viscosity
Dv = eta*h^3/3;                  % Viscous flexural rigidity (Pa m^3 s)

%% Profile specifications
dx = 50;
L = 40e3; % length of ice shelf 
x = 0:dx:L; 
x_km = x/1000;


%% solving viscous profile

% calculate beta 
beta = (rho_w*g/(omega*Dv))^(1/4);

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
