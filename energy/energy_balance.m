%% Test location in Thwaites region
lat = -75.5;
lon = -106.5;
altitude = 18e3;

%% Energy Generation

date_start = datenum(2026, 1, 18, 0, 0, 0);

% Simulation duration
n_days = 14;

% 5-minute resolution
dt_minutes = 5;
time_num = date_start : dt_minutes/(24*60) : date_start + n_days;

% Time in hours
time_hours = (time_num - date_start)*24;

% Preallocate
P_solar = zeros(size(time_num));

% Calculate solar power
for i = 1:length(time_num)

    [P_solar(i), ~, ~, ~, ~] = ...
        calculate_solar_power(lat, lon, time_num(i));

end

% Integrate power to get cumulative energy

% P_solar is W, time_hours is hours
% W*hr / 1000 = kWh
E_solar = cumtrapz(time_hours, P_solar) / 1000;

% Total solar energy over the day
E_total = E_solar(end);

fprintf('Total solar energy = %.2f kWh\n', E_total);


% Plot

figure;

yyaxis left
plot(time_hours, P_solar, 'LineWidth', 2);
ylabel('Solar Power (W)');

yyaxis right
plot(time_hours, E_solar, 'LineWidth', 2);
ylabel('Cumulative Solar Energy (kWh)');

xlabel('Time of Day UTC (hours)');
title(sprintf('Solar Power and Cumulative Solar Energy — %s', ...
    datestr(date_start, 'mmmm dd, yyyy')));


grid on;
box on;


%% Plot power usage profile
peak_draw = 2054.6; % W

V_ground = 20; % m/s
V_wind = 0; % m/s


[P_el, V_air, Cl, Cd, Ft] = calculate_flight_power(V_ground, V_wind);

% Constant electrical power usage
P_usage = P_el * ones(size(time_num));

% Cumulative energy usage [kWh]
E_usage = cumtrapz(time_hours, P_usage) / 1000;

% Total daily energy usage
fprintf('Aircraft power usage = %.2f W\n', P_el);
fprintf('Total daily energy usage = %.2f kWh\n', E_usage(end));

% Plot
figure;

yyaxis left
plot(time_hours, P_usage, 'LineWidth', 2);
ylabel('Power Usage (W)');
yline(peak_draw, 'r--', sprintf('Peak Power Draw: %.1f W', peak_draw), 'LabelHorizontalAlignment','left');

yyaxis right
plot(time_hours, E_usage, 'LineWidth', 2);
ylabel('Cumulative Energy Usage (kWh)');

xlabel('Time of Day (hours)');
title('Aircraft Power and Cumulative Energy Usage');

grid on;
box on;


%% Net Power

battery_capacity = 27.2; % kWh

% Net power at each time
P_net = P_solar - P_el;

% Preallocate battery energy
E_battery = zeros(size(time_hours));

% Start fully charged
E_battery(1) = battery_capacity;

for i = 2:length(time_hours)

    % Time step [hr]
    dt = time_hours(i) - time_hours(i-1);

    % Net energy added during this step [kWh]
    dE = 0.5 * (P_net(i) + P_net(i-1)) * dt / 1000;

    % Update battery from previous timestep
    E_battery(i) = E_battery(i-1) + dE;

    % Enforce storage limits immediately
    E_battery(i) = min(E_battery(i), battery_capacity);
    E_battery(i) = max(E_battery(i), 0);

end

fprintf('Net energy over the day = %.2f kWh\n', E_battery(end));

% Plot

figure;

yyaxis left

plot(time_hours / 24, P_net, 'LineWidth', 2);
hold on;

% Zero-power line
yline(0, '--', 'LineWidth', 1);

ylabel('Net Power (W)');

yyaxis right

plot(time_hours /24, E_battery, 'LineWidth', 2);
ylabel('Cumulative Net Energy (kWh)');
yline(battery_capacity, 'k-', ...
    sprintf('Battery Capacity (%.1f kWh)', battery_capacity), ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Battery Capacity');
xlabel('Time (Days)');

title(sprintf('Net Power and Cumulative Energy — %s to %s', ...
    datestr(date_start, 'mmmm dd, yyyy'), datestr(date_start + n_days, 'mmmm dd, yyyy')));

% xlim([0 24]);
% xticks(0:2:24);

grid on;
box on;
% 
% exportgraphics(gcf, ...
%     fullfile(figure_dir, 'two_week.png'), ...
%     'Resolution', 300);

%% Cleaned up plot

figure; 

% --- Power ---
subplot(2,1,1);
hold on;

plot(time_hours, P_solar, 'b:',  'LineWidth', 2);
plot(time_hours, -P_usage, 'b--', 'LineWidth', 2);
plot(time_hours, P_net,   'b-', 'LineWidth', 2);

yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
title('Instantaneous Power')
ylabel('Power (W)');
legend('Power Generated', 'Power Used', 'Net Power', ...
       'Location','southeast');
grid on;


% --- Battery energy ---
subplot(2,1,2);
battery_percent = 100 * E_battery / battery_capacity;
plot(time_hours, battery_percent, ...
    'Color', [0.8500 0.3250 0.0980], ...
    'LineWidth', 2);

title('Battery Capacity Percentage')

ylabel('Percentage of Battery Capacity');
xlabel('Time (hours)');
ylim([0 110]);
grid on;

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/EnergyBudget';

% Create folder if it doesn't exist
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

exportgraphics(gcf, ...
    fullfile(figure_dir, 'single_day.png'), ...
    'Resolution', 300);


%% plot varying air speeds
peak_draw = 2054.6; % W

V_ground = 5:5:40; % m/s
V_wind = 0; % m/s

figure;
hold on;

% Transparency range
alpha_min = 0.2;
alpha_max = 1.0;

for j = 1:length(V_ground)

    v = V_ground(j);

    [P_el, V_air, Cl, Cd, Ft] = ...
        calculate_flight_power(v, V_wind);

    % Constant electrical power usage
    P_usage = P_el * ones(size(time_num));

    % Cumulative energy usage [kWh]
    E_usage = cumtrapz(time_hours, P_usage) / 1000;

    % Increasing transparency with ground speed
    alpha = alpha_min + ...
        (alpha_max - alpha_min)*(j-1)/(length(V_ground)-1);

    % Power
    yyaxis left

    plot(time_hours, P_usage, '-', ...
        'Color', [0 0.4470 0.7410 alpha], ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('%d m/s', v));

    ylabel('Power Usage (W)');

    % Energy
    yyaxis right

    plot(time_hours, E_usage, '-', ...
        'Color', [0.8500 0.3250 0.0980 alpha], ...
        'LineWidth', 2, ...
        'HandleVisibility', 'Off');

    ylabel('Cumulative Energy Usage (kWh)');

end

% Peak draw line
yyaxis left

yline(peak_draw, 'r--', ...
    sprintf('Peak Power Draw: %.1f W', peak_draw), ...
    'LabelHorizontalAlignment', 'left', ...
    'LineWidth', 1.5, ...
    'HandleVisibility', 'Off');

% Formatting

xlabel('Time of Day (hours)');
title('Aircraft Power and Cumulative Energy Usage');

xlim([0 24]);
xticks(0:2:24);

legend('Location','best');

grid on;
box on;

exportgraphics(gcf, ...
    fullfile(figure_dir, 'varying_airspeed.png'), ...
    'Resolution', 300);

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