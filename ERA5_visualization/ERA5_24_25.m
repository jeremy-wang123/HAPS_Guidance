%% Playing around with ERA5 dataset

%{
Climate Reanalysis Data
- Area: Latrange: -60:-90, Lonrange: -180:180
- Times: 24 hours
- Years: 2024-2025
- Months: Dec, Jan, Feb
- Variables: u and v wind components
- Pressure: 70hPA ~ 18km
%} 

addpath('/Users/jeremywang/Documents/MATLAB/ERA5');

filename = 'data_stream-oper_stepType-instant.nc';

lon = ncread(filename,'longitude');
lat = ncread(filename,'latitude');
pressure = ncread(filename,'pressure_level');
time_raw = ncread(filename,'valid_time');

% remove pressure column
u = squeeze(ncread(filename, 'u'));
v = squeeze(ncread(filename, 'v'));

% Convert Unix time to MATLAB datetime
time = datetime(1970,1,1,0,0,0) + seconds(time_raw);

% Wind speed
wind_speed = hypot(u,v);

%% Coordinate conversion
% Full ERA5 coordinate grid
[Lon,Lat] = meshgrid(lon,lat);

% Convert ERA5 coordinates to Antarctic polar stereographic
[x_wind,y_wind] = ll2ps(Lat,Lon);

% Thwaites map limits (polar stereographic)
mapzoomps('Thwaites Glacier');
xl = xlim;
yl = ylim;

x_min = xl(1);
x_max = xl(2);

y_min = yl(1);
y_max = yl(2);
close

% Select only ERA5 points inside Thwaites map extent
inside = ...
    x_wind >= x_min & x_wind <= x_max & ...
    y_wind >= y_min & y_wind <= y_max;

% Determine local east/north directions
d = 1e-4;

[x_e,y_e] = ll2ps(Lat,Lon + d);

ex = x_e - x_wind;
ey = y_e - y_wind;

e_mag = hypot(ex,ey);

ex = ex ./ e_mag;
ey = ey ./ e_mag;


[x_n,y_n] = ll2ps(Lat + d,Lon);

nx = x_n - x_wind;
ny = y_n - y_wind;

n_mag = hypot(nx,ny);

nx = nx ./ n_mag;
ny = ny ./ n_mag;


%% select time step to plot

it = 1;

% wind vectors for selected time value
u_plot = u(:,:,it);
v_plot = v(:,:,it);

% ERA5 arrays are lon x lat, so transpose
u_map = u_plot';
v_map = v_plot';

% Convert geographic winds to projected x/y components
u_ps = u_map .* ex + v_map .* nx;
v_ps = u_map .* ey + v_map .* ny;

%
figure;

measuresps('speed');
mapzoomps('Thwaites Glacier')
hold on

quiver( ...
    x_wind(inside), ...
    y_wind(inside), ...
    u_ps(inside), ...
    v_ps(inside), ...
    'k', ...
    'LineWidth',1.2, ...
    'AutoScale','on', ...
    'AutoScaleFactor',1);

title(sprintf('ERA5 Winds over Thwaites — %s', ...
    string(time(it))));




%% Animation

% it_array = 1:24; % Jan 1, 24 hours
% it_array = 1:24:(1 + 24*29); % 2024, Jan
% 
target_hour = 12;
it_array = find( ...
    time >= datetime(2024,12,1) & ...
    time <  datetime(2025,3,1) & ...
    hour(time) == target_hour); % Dec 2024 - Jan 2025

% quiver_scale
quiver_scale = 1000;

% Create map once
figure;

measuresps('speed');
mapzoomps('Thwaites Glacier');
hold on;

% First time step
it = it_array(1);

u_map = u(:,:,it)';
v_map = v(:,:,it)';

u_ps = u_map .* ex + v_map .* nx;
v_ps = u_map .* ey + v_map .* ny;

% Create quiver object
q = quiver( ...
    x_wind(inside), ...
    y_wind(inside), ...
    quiver_scale*u_ps(inside), ...
    quiver_scale*v_ps(inside), ...
    'k', ...
    'LineWidth', 1.2, ...
    'AutoScale', 'on', ...
    'AutoScaleFactor', 1);

% Reference vector with white box — top right
xl = xlim;
yl = ylim;

V_ref = 20;   % m/s


% Box size
box_width  = 0.18 * range(xl);
box_height = 0.1 * range(yl);

% Box position — top right with margin
box_x = xl(2) - box_width  - 0.03*range(xl);
box_y = yl(2) - box_height - 0.03*range(yl);

% White background box
rectangle( ...
    'Position',[box_x, box_y, box_width, box_height], ...
    'FaceColor','w', ...
    'EdgeColor','k', ...
    'LineWidth',1);

% Center of box
box_x_center = box_x + box_width/2;
box_y_center = box_y + box_height/4;

% Arrow length
arrow_length = quiver_scale * V_ref;

% Start arrow so its midpoint is at box center
x_ref = box_x_center - arrow_length/2;
y_ref = box_y_center;

quiver( ...
    x_ref, y_ref, ...
    arrow_length, 0, ...
    0, ...
    'k', ...
    'LineWidth',1.5, ...
    'MaxHeadSize',0.8);

% Label above arrow
text( ...
    box_x + 0.5*box_width, ...
    box_y + 0.75*box_height, ...
    sprintf('%g m/s',V_ref), ...
    'FontSize',8, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle');

title(sprintf('ERA5 Winds over Thwaites — %s', ...
    string(time(it))));

% Animate remaining timesteps
for i = 2:length(it_array)

    it = it_array(i);

    % Extract wind
    u_map = u(:,:,it)';
    v_map = v(:,:,it)';

    % Rotate into polar stereographic coordinates
    u_ps = u_map .* ex + v_map .* nx;
    v_ps = u_map .* ey + v_map .* ny;

    % Update existing arrows
    q.UData = quiver_scale*u_ps(inside);
    q.VData = quiver_scale*v_ps(inside);
    
    % Reference vector with white box — top right
    xl = xlim;
    yl = ylim;
    
    V_ref = 20;   % m/s
    
    
    % Box size
    box_width  = 0.18 * range(xl);
    box_height = 0.1 * range(yl);
    
    % Box position — top right with margin
    box_x = xl(2) - box_width  - 0.03*range(xl);
    box_y = yl(2) - box_height - 0.03*range(yl);
    
    % White background box
    rectangle( ...
        'Position',[box_x, box_y, box_width, box_height], ...
        'FaceColor','w', ...
        'EdgeColor','k', ...
        'LineWidth',1);
    
    % Center of box
    box_x_center = box_x + box_width/2;
    box_y_center = box_y + box_height/4;
    
    % Arrow length
    arrow_length = quiver_scale * V_ref;
    
    % Start arrow so its midpoint is at box center
    x_ref = box_x_center - arrow_length/2;
    y_ref = box_y_center;
    
    quiver( ...
        x_ref, y_ref, ...
        arrow_length, 0, ...
        0, ...
        'k', ...
        'LineWidth',1.5, ...
        'MaxHeadSize',0.8);
    
    % Label above arrow
    text( ...
        box_x + 0.5*box_width, ...
        box_y + 0.75*box_height, ...
        sprintf('%g m/s',V_ref), ...
        'FontSize',8, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');

    % Update title
    title(sprintf('ERA5 Winds over Thwaites — %s', ...
        string(time(it))));

    % Update figure
    drawnow;

    % Controls animation speed
    pause(0.5);
end

%% Save animation
output_folder = ...
    ['/Users/jeremywang/Library/CloudStorage/' ...
     'GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/' ...
     'ERA5_frames/Daily'];

if ~exist(output_folder,'dir')
    mkdir(output_folder);
end

% Create figure once
fig = figure('Position',[100 100 900 750]);

measuresps('speed');
mapzoomps('Thwaites Glacier');
hold on;

ax = gca;

% Initialize quiver
it = it_array(1);

u_map = u(:,:,it)';
v_map = v(:,:,it)';

u_ps = u_map .* ex + v_map .* nx;
v_ps = u_map .* ey + v_map .* ny;

quiver_scale = 1000;

q = quiver( ...
    x_wind(inside), ...
    y_wind(inside), ...
    quiver_scale*u_ps(inside), ...
    quiver_scale*v_ps(inside), ...
    0, ...
    'k', ...
    'LineWidth',1);

% Generate frames
for i = 1:length(it_array)

    it = it_array(i);

    u_map = u(:,:,it)';
    v_map = v(:,:,it)';

    u_ps = u_map .* ex + v_map .* nx;
    v_ps = u_map .* ey + v_map .* ny;

    q.UData = quiver_scale*u_ps(inside);
    q.VData = quiver_scale*v_ps(inside);
    % Reference vector with white box — top right
    xl = xlim;
    yl = ylim;
    
    V_ref = 20;   % m/s
    
    
    % Box size
    box_width  = 0.18 * range(xl);
    box_height = 0.1 * range(yl);
    
    % Box position — top right with margin
    box_x = xl(2) - box_width  - 0.03*range(xl);
    box_y = yl(2) - box_height - 0.03*range(yl);
    
    % White background box
    rectangle( ...
        'Position',[box_x, box_y, box_width, box_height], ...
        'FaceColor','w', ...
        'EdgeColor','k', ...
        'LineWidth',1);
    
    % Center of box
    box_x_center = box_x + box_width/2;
    box_y_center = box_y + box_height/4;
    
    % Arrow length
    arrow_length = quiver_scale * V_ref;
    
    % Start arrow so its midpoint is at box center
    x_ref = box_x_center - arrow_length/2;
    y_ref = box_y_center;
    
    quiver( ...
        x_ref, y_ref, ...
        arrow_length, 0, ...
        0, ...
        'k', ...
        'LineWidth',1.5, ...
        'MaxHeadSize',0.8);
    
    % Label above arrow
    text( ...
        box_x + 0.5*box_width, ...
        box_y + 0.75*box_height, ...
        sprintf('%g m/s',V_ref), ...
        'FontSize',8, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
   

    title(sprintf('ERA5 Winds over Thwaites — %s', ...
        string(time(it))));

    drawnow;

    filename = fullfile( ...
        output_folder, ...
        sprintf('frame_%03d.png',i));

    print(fig,filename,'-dpng','-r150');

end
%% Save as a video
frame_folder = ...
    ['/Users/jeremywang/Library/CloudStorage/' ...
     'GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/' ...
     'ERA5_frames/Daily'];

% Output MP4
video_file = ...
    ['/Users/jeremywang/Library/CloudStorage/' ...
     'GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/' ...
     'ERA5/ERA5_Thwaites_Daily.mp4'];

video = VideoWriter(video_file,'MPEG-4');

video.FrameRate = 2;
video.Quality = 100;

open(video);

% Find all frames
files = dir(fullfile(frame_folder,'frame_*.png'));

% Sort by filename
[~,idx] = sort({files.name});
files = files(idx);

% Add PNGs to video
for i = 1:length(files)

    filename = fullfile(frame_folder,files(i).name);

    img = imread(filename);

    writeVideo(video,img);

end

close(video);

fprintf('Video saved to:\n%s\n',video_file);
