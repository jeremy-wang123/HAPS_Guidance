%% Test location in Thwaites region

lat = -75.5;
lon = -106.5;

[x, y] = ll2ps(lat, lon);

fprintf('x = %.1f m\n', x);
fprintf('y = %.1f m\n', y);

% Time
time_num = datenum(2026, 1, 15, 12, 0, 0);
time_zone = 0;

% Coordinate-system rotation relative to north
rotation = 0;

% No daylight savings time
dst = false;

% Solar position

[angles, projection] = solarPosition( ...
    time_num, ...
    lat, ...
    lon, ...
    time_zone, ...
    rotation, ...
    dst);

% Outputs
zeta = angles(:,1);       % solar zenith angle [deg] 
solar_az = angles(:,2);   % solar azimuth [deg]

fprintf('Solar zenith angle = %.2f deg\n', zeta);
fprintf('Solar azimuth      = %.2f deg\n', solar_az);

% Level aircraft assumption

theta_incidence = zeta;

solar_projection = max(cosd(theta_incidence),0);

fprintf('cos(incidence angle) = %.3f\n', solar_projection);