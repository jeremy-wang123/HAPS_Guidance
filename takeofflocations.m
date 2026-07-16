%% takeoff locations

% visualize takeoff locations at Union Glacier station, McMurdo Station,
% and Rothera Research Station

% Station coordinates [degrees]
station_names = {'Union Glacier', 'McMurdo Station', 'Rothera Station'};

lat = [...
    -79.7680;   % Union Glacier Camp
    -77.8463;   % McMurdo Station
    -67.568889];  % Rothera Station

lon = [...
    -83.2617;   % Union Glacier Camp
     166.6682;  % McMurdo Station
    -68.1248];  % Rothera Station

figure;
% Create a full-continent Antarctic map
mapzoomps('Rothera Station');
xl = xlim;
yl = ylim;
measuresps('speed')
hold on;

% Antarctic grounding line/coastline
measuresps('gl', 'k');

% Convert station coordinates to polar stereographic x-y coordinates
[x_station, y_station] = ll2ps(lat, lon);

% Plot station markers
scatter(x_station, y_station, 100, ...
    'filled', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2);

% Add station labels
text(x_station, y_station, station_names, ...
    'FontSize', 11, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left');

title('Antarctic Stations');

axis image;
xlim(xl);
ylim(yl);


%% Concentric plots around takeoff points

speed = 20; % m/s (estimation)
v = speed*(60*60*24); % m/day

takeoff_loc = [-79.7680 -83.2617]; % Union Glacier Camp
[x_station, y_station] = ll2ps(takeoff_loc(1), takeoff_loc(2));

campaigns = [0.25 0.5 1 1.5 2 3 4 5];

figure;

measuresps('speed', 'log');
colormap(parula);
cb = colorbar;
set(gca,'ColorScale');
cb.Ticks = [0 1 2 3];
cb.TickLabels = {'1','10','100','1000'};
cb.Label.String = 'Speed (m/yr)';
hold on;

% plotting point
scatter(x_station, y_station, 100, ...
    "red",...
    'filled', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2, ...
    'DisplayName','Union Glacier');

% plotting 1-7 day campaign
for i=campaigns
    campaign = i; 
    max_distance = v*campaign/2; % m

    % Generate concentric circles around the takeoff location
    theta = linspace(0, 2*pi, 100); % Angle for circle
    x_circle = x_station + (max_distance * cos(theta)); % Longitude
    y_circle = y_station + (max_distance * sin(theta)); % Latitude

    plot(x_circle, y_circle, ...
    'DisplayName', sprintf('%g-Day Campaign', i));
end
legend('show', 'Location', 'southeast');
axis equal
xlabel('x (m)');
ylabel('y (m)');
title('Max Distance from Union Glacier per Campaign')

hold off;

%% Determine Swath Dimensions

beta_el = 40; % elevation beam width
beta_az = 40; % azimuth beam width
theta = 40; % look angle
z = 20; % altitude of flight

S_range = z * (tand(theta + beta_el/2) - tand(theta - beta_el/2));
S_azimuth = (2*z / cosd(theta)) * tand(beta_az/2);

S_range_approx = (z*deg2rad(beta_el)/(cosd(theta)^2));
S_azimuth_approx = (z*deg2rad(beta_az)/cosd(theta));

% instantaneous area
S_area = pi/4 * S_range * S_azimuth;
S_area_approx = pi/4 * S_range_approx * S_azimuth_approx;

%% NISAR Cross checking calculations
beta_el = 12; % elevation beam width
beta_az = 40; % azimuth beam width
theta = 37; % look angle
z = 747; % altitude of flight

S_range = z * (tand(theta + beta_el/2) - tand(theta - beta_el/2));
S_range_approx = (z*deg2rad(beta_el)/(cosd(theta)^2));
