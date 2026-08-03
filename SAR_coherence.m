%% Calculations for the uncertainty and coherence of radar

% constants
c = 299792458; % speed of light, m/s
f = 1.2e9; % cycles/sec
theta = 40; % look angle in degrees
lambda = c/f; % meters

r_slant = 0.75; % m from Gamma manual
beta_az = 40; % azimuth beam width
altitude = 18e3; % m

ry = r_slant / sind(theta); % range resolution (m)
rx = 0.5; % azimuth resolution (m) (listed between 0.2 - 0.5 m)

%% Exponential Temporal coherence
%{
Assumptions made within calculation: 
1) Temporal coherence, sy, and sz are functions of time. Temporal coherence
follows a first order decaying exponential, where 5*tau occurs at 6
days, when the temporal coherence is approximately 0. (gamma(6 days) =
0.7)
2) Assume linear relationship between sy and sz, where sz = k*sy
%}


% parameters
tau = 6/5 * 24 * 3600; % seconds
k = [0.25 1 4]; % linear coefficient relating sy and sz
t = 0:1:(3600*24*10); % 10 day simulation period
t_days = t/(24*3600);

% plot gamma as a function of time
figure;
% Calculate the temporal coherence function
gamma = exp(-t/tau);
plot(t_days,gamma);
xlabel('Time (days)');
ylabel('\gamma_{temporal}');
title('Decaying Exponential \gamma_{temporal}')

% plotting sy and sz
figure;
hold on;
grid on;
box on;

colors = lines(length(k));

for i = 1:length(k)

    sy = (lambda/(4*pi)) .* sqrt(2*t./tau) ./ ...
         sqrt(sind(theta).^2 + k(i)^2 .* cosd(theta).^2);

    sz = k(i) .* sy;

    % Horizontal displacement (solid)
    plot(t_days, sy, '-', ...
        'Color', colors(i,:), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('s_y, k = %.2f', k(i)));

    % Vertical displacement (dashed)
    plot(t_days, sz, '--', ...
        'Color', colors(i,:), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('s_z, k = %.2f', k(i)));

end

xlabel('Time (days)');
ylabel('RMS Scatterer Displacement (m)');
title('Estimated Scatterer Displacement vs Time');

legend('Location','northwest');

%% Linear Temporal Coherence
% Parameters for linear temporal coherence

T_zero = 6 * 24 * 3600;          % coherence reaches zero at 6 days, s
k = [0.25 1 4];                  % sz = k*sy

t = 0:60:(10 * 24 * 3600);       % 10-day simulation, 1-minute spacing
t_days = t/(24*3600);

% Linear temporal coherence
gamma_temporal = max(1 - t/T_zero, 0);

figure;
plot(t_days, gamma_temporal, 'LineWidth', 2);
grid on;
box on;

xlabel('Time (days)');
ylabel('\gamma_{temporal}');
title('Linear Temporal Coherence Decay');

xlim([0 10]);
ylim([0 1]);
xline(6, '--', 'Complete decorrelation');

% Estimate sy and sz from temporal coherence

figure;
hold on;
grid on;
box on;

colors = lines(length(k));

% The inverse equation is undefined at gamma = 0
valid = gamma_temporal > 0;

for i = 1:length(k)

    sy = NaN(size(t));
    sz = NaN(size(t));

    sy(valid) = (lambda/(4*pi)) .* ...
        sqrt( ...
        -2 .* log(gamma_temporal(valid)) ./ ...
        (sind(theta).^2 + k(i).^2 .* cosd(theta).^2) ...
        );

    sz(valid) = k(i) .* sy(valid);

    % Horizontal displacement spread
    plot(t_days, sy, '-', ...
        'Color', colors(i,:), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('s_y, k = %.2f', k(i)));

    % Vertical displacement spread
    plot(t_days, sz, '--', ...
        'Color', colors(i,:), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('s_z, k = %.2f', k(i)));

end

xline(6, 'k:', 'Complete decorrelation', ...
    'HandleVisibility', 'off');

xlabel('Time (days)');
ylabel('RMS Scatterer Displacement (m)');
title('Scatterer Displacement for Linear Coherence Decay');

xlim([0 6]);
legend('Location', 'northwest');

%% Temporal coherence vs horizontal RMS scatterer displacement
% Assumption: sz = k*sy

k = [0.25 1 4];

sy = linspace(0,0.12,500);    % horizontal RMS displacement, m
colors = lines(length(k));

figure;
hold on;
grid on;
box on;

for i = 1:length(k)

    sz = k(i).*sy;

    gamma_temporal = exp( ...
        -0.5 .* (4*pi/lambda).^2 .* ...
        (sy.^2 .* sind(theta).^2 + ...
         sz.^2 .* cosd(theta).^2) );

    plot(sy, gamma_temporal, ...
        'LineWidth',2, ...
        'Color',colors(i,:), ...
        'DisplayName',sprintf('k = %.2f',k(i)));
end

xlabel('Horizontal RMS Scatterer Displacement, s_y (m)');
ylabel('\gamma_{temporal}');
title('Temporal Coherence vs. RMS Scatterer Displacement');

xlim([0 max(sy)]);
ylim([0 1]);

legend('Location','northeast');

gamma_temporal = 0.9; % guess
%% Thermal coherence
% we can assume that gamma ~ 1 for SNR >> 1
% -15 dB for noise equivalent sigma zero
% estiamte something for sigma naught

SNR = 30; % complete guess, can't find the value
gamma_thermal = SNR / (SNR + 1);

%% spatial coherence
b_perp = 100; % perpendicular component to baseline vector (m)
R = altitude / cosd(theta); % slant range (m)

Bc = (lambda * R) / (2*ry * cosd(theta));
gamma_spatial = 1 - (2*b_perp*ry*cosd(theta) / (lambda*R));

%% rotational coherence
dphi_deg = 5; % difference in heading directions
dphi = deg2rad(dphi_deg);

dphi_c = rad2deg(lambda / (2*rx * sind(theta))); % maximum heading difference allowed
gamma_rotation = 1 - (2*rx*dphi*sind(theta))/lambda; % rotational coherence


%% total coherence

gamma = gamma_thermal * gamma_spatial * gamma_rotation * gamma_temporal; % total coherence

% L = 1 can be the lower bound, but need to determine the number of pixels
% system averages over (based on processing setup
L = 100; % number of looks (need to figure out)

sigma_phi = (sqrt(1-gamma^2) / (gamma * sqrt(2*L)));
sigma_z = (lambda/(4*pi*cosd(theta))) * sigma_phi;

