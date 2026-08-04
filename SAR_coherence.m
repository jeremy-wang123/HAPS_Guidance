%% Calculations for the uncertainty and coherence of radar

% constants
c = 299792458; % speed of light, m/s
f = 1.2e9; % cycles/sec
lambda = c/f; % meters


%% Linear Temporal Coherence
% Parameters for linear temporal coherence

T_zero = 6 * 24 * 3600;          % coherence reaches zero at 6 days, s
k = 1;                  % sz = k*sy

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


% The inverse equation is undefined at gamma = 0
valid = gamma_temporal > 0;

sy = NaN(size(t));
sz = NaN(size(t));

sy(valid) = (lambda/(4*pi)) .* ...
    sqrt( ...
    -2 .* log(gamma_temporal(valid)) ./ ...
    (sind(theta).^2 + k.^2 .* cosd(theta).^2) ...
    );

sz(valid) = k .* sy(valid);

% Horizontal displacement spread
plot(t_days, sy, '-', ...
    'LineWidth', 2);

xline(6, 'k:', 'Complete decorrelation', ...
    'HandleVisibility', 'off');

xlabel('Time (days)');
ylabel('RMS Scatterer Displacement (m)');
title('Scatterer Displacement for Linear Coherence Decay');

xlim([0 6]);

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
%% function for calculating sigma_z
L = 100; 
lambda = c/f;
dphi = 1; 
NESZ_dB = -20; 
sigma0_dB = -10;
b_perp = 100; 
time = 12*3600; % seconds

estimate_insar_uncertainty(L, lambda, dphi, NESZ_dB, sigma0_dB, b_perp, time)