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
