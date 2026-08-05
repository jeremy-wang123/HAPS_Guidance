%% Analytical solution to elastic beam bending

%% Parameter values 
%{
w0 = tidal amplitude at the edge of the ice shelf (determined from the
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

% extract key params
[h_shelf_mean,h_shelf_eff,h_front_mean,h_front_eff] = extract_shelf_thickness('Thwaites Glacier');
mean_amp = calculate_mean_amp('Thwaites', 'top2');

% parameters 
rho = 1025; % kg/m^3
g = 9.81; % m/s^2
h = h_front_eff; % use the mean effective thickness of calving front
w0 = mean_amp; % m
sigma_z = 0.002;  % m (vertical displacement uncertainty)

% parameters for elastic beam bending
I = h^3/12;
mu = 0.325; % poisson's ratio
E = 9.33e9; % Pa
D = (E*I)/(1-mu^2); % flexural rigidity

% parameters for viscous beam bending
T_hours = 24;            % Tidal period (hours), M2 example
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

%% Solving w(x) for elastic profile

lambda = (rho*g*(1-mu^2)/(4*E*I))^(1/4);
min_L = 5*pi/(4*lambda); % minimum ice shelf length for elastic estimation to make sense (according to Holdsworth)
w = w0 * (1 - (exp(-lambda*x) .* (cos(lambda*x) + sin(lambda*x))));
%% solving viscous profile

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
%% Settling-distance calculations

% calculate settling distance analytically
elastic_settling_distance = calculate_settling_distance(w, x, w0, sigma_z);

% calculate settling distance
viscous_settling_distance = calculate_settling_distance(w_peak, x, w0, sigma_z);


%% Plotting Elastic and Viscous Profiles

figure;
hold on;
plot(x_km, w, 'LineWidth', 2); % elastic
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
xline(elastic_settling_distance, ':', ...
    sprintf('Elastic Settling distance = %.1f km', elastic_settling_distance), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left');

xline(viscous_settling_distance, ':', ...
    sprintf('Viscous Settling distance = %.1f km', viscous_settling_distance), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Deflection at Peak Tide Thwaites Glacier');

xlim([0, L/1000]);

% Include the displacement profile and uncertainty bounds in the y-limits
all_y_values = [w(:); upper_bound; lower_bound];

w_range = max(all_y_values) - min(all_y_values);

% Avoid zero padding if the plotted values are nearly constant
if w_range == 0
    vertical_padding = max(abs(all_y_values))*0.1;
else
    vertical_padding = 0.20*w_range;
end

ylim([min(all_y_values) - vertical_padding, ...
      max(all_y_values) + vertical_padding]);

legend('Elastic', 'Viscous', 'Location','southeast')
grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);

% figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending';
% exportgraphics(gcf, fullfile(figure_dir,'elastic_viscous_beam_bending_thwaites.jpg'), ...
% 'Resolution',300);

%% plotting curvature at peak tide

% analytically solved elastic curvature
elastic_curvature = 2*w0*lambda^2 * exp(-lambda*x) .* (cos(lambda*x) - sin(lambda*x));

% viscous curvature
viscous_slope = gradient(w_peak, x);          % dw/dx
viscous_curvature = gradient(viscous_slope, x);  % d^2w/dx^2

figure;
hold on;
plot(x_km, viscous_curvature, 'LineWidth', 2);
plot(x_km,elastic_curvature, 'LineWidth', 2);

xlabel('Distance from grounding line (km)');
ylabel('Curvature, $d^2w/dx^2 (m^{-1})$', 'Interpreter', 'latex');
title('Curvature at peak tide');

yline(0, 'k--');
xlim([0, L/1000]);

legend('elastic','viscous')
grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);


% figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending';
% exportgraphics(gcf, fullfile(figure_dir,'beam_curvature.jpg'), ...
% 'Resolution',300);

%% Deflection at lf

% flexural wave-length (features smaller than lf are support by rigidity,
% features larger than lf are at isostatic equilibrium)
lf = (4*D/(rho*g))^(1/4);

% Find closest grid point
[~,idx] = min(abs(x-lf));

% Time vector (two tidal cycles)
Nt = 500;
t = linspace(0,2*T,Nt);

% Initialize
w_lf = zeros(size(t));

% Displacement at x = lf
for i = 1:Nt

    % Instantaneous displacement profile
    w = imag(W .* exp(1i*omega*t(i)));

    % Record displacement at x = lf
    w_lf(i) = w(idx);

end

% Applied tidal forcing
w_tide = w0*sin(omega*t);

% Plot
figure;
hold on;

plot(t/3600, w_lf, 'LineWidth', 2, ...
    'DisplayName','viscous');

plot(t/3600, w_tide, '--k', 'LineWidth', 2, ...
    'DisplayName','$w_{tide}$');

xlabel('Time (hours)');
ylabel('Vertical displacement (m)');
title(sprintf('Response at $x=\\ell_f = %.2f$ km', lf/1000), ...
    'Interpreter','latex');

legend('Interpreter','latex','Location','best');

grid on;
box on;
set(gca,'FontSize',12,'LineWidth',1);


% figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending';
% exportgraphics(gcf, fullfile(figure_dir,'viscous_phase_lag.jpg'), ...
% 'Resolution',300);