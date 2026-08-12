%% Preliminary plots

%% Quiver and color plot of ice velocities
figure
measuresps('speed','log','alpha',0.7)
cb = colorbar;
cb.Ticks = [0 1 2 3];
cb.TickLabels = {'1','10','100','1000'};
hold on
measuresps('gl','r')

% manually plotting quivers because measuresps doesn't show all velocities
% Get current map limits
xl = xlim;
yl = ylim;

% Read MEaSUREs velocity data
x = ncread('antarctica_ice_velocity_450m_v2.nc','x');
y = ncread('antarctica_ice_velocity_450m_v2.nc','y');
vx = ncread('antarctica_ice_velocity_450m_v2.nc','VX');
vy = ncread('antarctica_ice_velocity_450m_v2.nc','VY');

% Downsample manually 
skip = 150;   

xq = x(1:skip:end);
yq = y(1:skip:end);
vxq = vx(1:skip:end,1:skip:end);
vyq = vy(1:skip:end,1:skip:end);

[X,Y] = meshgrid(xq,yq);

h = quiver(X,Y,vxq',vyq','k','LineWidth',0.8);
xlabel('x (m)');
ylabel('y (m)');
xlim(xl)
ylim(yl)
set(h,'AutoScaleFactor',1)
axis equal

%% Quiver plot for Thwaites
figure
mapzoomps('Pine Island Glacier')
measuresps('speed','log','alpha',0.7)
cb = colorbar;
cb.Ticks = [0 1 2 3];
cb.TickLabels = {'1','10','100','1000'};
hold on
measuresps('gl','r')

% manually plotting quivers because measuresps doesn't show all velocities
% Get current map limits
xl = xlim;
yl = ylim;

% Read MEaSUREs velocity data
x = ncread('antarctica_ice_velocity_450m_v2.nc','x');
y = ncread('antarctica_ice_velocity_450m_v2.nc','y');
vx = ncread('antarctica_ice_velocity_450m_v2.nc','VX');
vy = ncread('antarctica_ice_velocity_450m_v2.nc','VY');

% Downsample manually 
skip = 30;   

xq = x(1:skip:end);
yq = y(1:skip:end);
vxq = vx(1:skip:end,1:skip:end);
vyq = vy(1:skip:end,1:skip:end);

[X,Y] = meshgrid(xq,yq);

h = quiver(X,Y,vxq',vyq','k','LineWidth',0.8);
title('Velocity plot of Thwaites Glacier')
xlim(xl)
ylim(yl)
xlabel('x (m)');
ylabel('y (m)');
set(h,'AutoScaleFactor',1)

%% Sinusoidal plot of tidal amplitude 14 day time series 
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008') 
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008'; 

% use lat and lon for the Thwaites Glacier
% We can switch to the polar coordinates later
latlim = [-76.0 -74.5];
lonlim = [-108 -103];

% specific coordinate choice 
lat = -75;
lon = -105; 

% time vector, hourly
t0 = datenum(2026,1,1,0,0,0);
time=t0+(0:(14*24))/24;

[z, conlist] = tmd_tide_pred(Model, time, lat, lon, 'z');
conList_clean = lower(strtrim(cellstr(conlist)));

figure;
plot(time - t0, z, 'b', 'LineWidth',1);
grid on;
xlabel('Time, days since 2026');
ylabel('Tidal height');
%% now plot each of the constit contributions
% specific coordinate choice 
lat = -75;
lon = -105; 

figure;

for i  =1:length(conList_clean)
    zi = tmd_tide_pred(Model, time, lat, lon, 'z', i);
    
    subplot(5,2,i)
    plot(time - t0, zi, 'LineWidth', 1, 'DisplayName', conList_clean{i})
    title(conList_clean{i})
    xlabel('Time, days since 2026');
    ylabel('Tidal height');
end


%% plot K1 and O1 sinusoidal contributions, averaged around Thwaites (from paper)
clear;clc;
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008') 
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008'; 

latlim = [-76.0 -74.5];
lonlim = [-108 -103];

lats = latlim(1):0.5:latlim(2);
lons = lonlim(1):0.5:lonlim(2);

% time vector, hourly
t0 = datenum(2026,1,1,0,0,0);
time=t0+(0:(46*24))/24;

% Check constituent list
[~,~,~,conList] = tmd_extract_HC(Model,lats(1),lons(1),'z');
conList_clean = lower(strtrim(cellstr(conList)));

idx_k1 = find(strcmp(conList_clean,'k1'));
idx_o1 = find(strcmp(conList_clean,'o1'));

k1_total = NaN(length(lats), length(lons), length(time));
o1_total = NaN(length(lats), length(lons), length(time));
z_total = NaN(length(lats), length(lons), length(time));
for i=1:length(lats)
    for j=1:length(lons)
        k1 = tmd_tide_pred(Model, time, lats(i), lons(j), 'z', idx_k1);
        o1 = tmd_tide_pred(Model, time, lats(i), lons(j), 'z', idx_o1);
        
        k1_total(i,j,:) = k1;
        o1_total(i,j,:) = o1;
        z_total(i,j,:) = k1 + o1;
    end
end


k1mean    = squeeze(mean(k1_total, [1 2], 'omitnan'));
o1mean    = squeeze(mean(o1_total, [1 2], 'omitnan'));
totalmean = squeeze(mean(z_total,  [1 2], 'omitnan'));

figure;
subplot(3,1,1)
plot(time - t0, k1mean)
title('K1 Thwaites Averaged')
xlabel('Days')
ylabel('Height')
subplot(3,1,2)
plot(time - t0, o1mean)
title('O1 Thwaites Averaged')
xlabel('Days')
ylabel('Height')
subplot(3,1,3)
plot(time - t0, totalmean)
title('K1+O1 Thwaites Averaged')
xlabel('Days')
ylabel('Height')
%% K1 + O1 amplitude map from CATS2008 near Thwaites

addpath('/Users/jeremywang/Documents/MATLAB/CATS2008')
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008';

% Thwaites region
latrange = -78.0:0.02:-73;
lonrange = -111:0.02:-99;

[lon,lat] = meshgrid(lonrange,latrange);

% Check constituent list
[~,~,~,conList] = tmd_extract_HC(Model,lat(1),lon(1),'z');
conList_clean = lower(strtrim(cellstr(conList)));

idx_k1 = find(strcmp(conList_clean,'k1'));
idx_o1 = find(strcmp(conList_clean,'o1'));

[KO_amp,~,~,~] = tmd_extract_HC(Model,lat,lon,'z',[idx_k1, idx_o1]);

KO_amp = squeeze(sum(KO_amp, 1));

% Convert CATS lon/lat grid to polar stereographic x/y used by CATS
[x,y] = ll2ps(lat,lon);


figure
% ice speed plot
ax1 = axes;
mapzoomps('Thwaites Glacier')
measuresps('speed','log')
measuresps('gl', 'k')
hold on

axis(ax1,'equal')
axis(ax1,'manual')

xl = xlim(ax1);
yl = ylim(ax1);

colormap(parula)
cb = colorbar('westoutside');
cb.Ticks = [0 1 2 3];
cb.TickLabels = {'1','10','100','1000'};
ylabel(cb,'Ice speed (m/yr)')

xlabel('x (m)');
ylabel('y (m)');
% tidal amplitude plot
ax2 = axes;

pcolor(ax2, x, y, KO_amp)
shading(ax2,'interp')
set(findobj(ax2,'Type','Surface'), ...
    'FaceAlpha',0.5, ...
    'EdgeColor','none')

colormap(ax2, turbo)
cb2 = colorbar(ax2,'eastoutside');
ylabel(cb2,'K1 + O1 amplitude (m)')

% Match axes exactly
set(ax2, ...
    'Position', ax1.Position, ...
    'Color','none', ...
    'XLim', xl, ...
    'YLim', yl, ...
    'DataAspectRatio', ax1.DataAspectRatio, ...
    'PlotBoxAspectRatio', ax1.PlotBoxAspectRatio, ...
    'XTick', [], ...
    'YTick', [], ...
    'Box','off')

linkaxes([ax1 ax2],'xy')

% Put ice-speed axes visually underneath
uistack(ax1,'bottom')

title(ax1,'CATS2008 K1 + O1 tidal amplitude over ice surface speed')

%% Now recreate figure 1; K1 and O1 for all Antarctica
clear;clc;

addpath('/Users/jeremywang/Documents/MATLAB/CATS2008')
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008';

latrange = -90.0:0.02:-60;
lonrange = -180:0.02:180;

[lon,lat] = meshgrid(lonrange,latrange);

% Check constituent list
[~,~,~,conList] = tmd_extract_HC(Model,lat(1),lon(1),'z');
conList_clean = lower(strtrim(cellstr(conList)));

idx_k1 = find(strcmp(conList_clean,'k1'));
idx_o1 = find(strcmp(conList_clean,'o1'));

% Extract harmonic amplitudes
[K1_amp,~,~,~] = tmd_extract_HC(Model,lat,lon,'z',idx_k1);
[O1_amp,~,~,~] = tmd_extract_HC(Model,lat,lon,'z',idx_o1);

% Combined diurnal amplitude
KO_amp = K1_amp + O1_amp;


% Convert CATS lon/lat grid to polar stereographic x/y used by CATS
[x,y] = ll2ps(lat,lon);


figure
% ice speed plot
ax1 = axes;
measuresps('speed', 'log')

hold on
measuresps('gl','k')

axis(ax1,'equal')
axis(ax1,'manual')

xl = xlim(ax1);
yl = ylim(ax1);

colormap(ax1, parula)
cb1 = colorbar(ax1,'westoutside');
cb1.Ticks = [0 1 2 3];
cb1.TickLabels = {'1','10','100','1000'};
ylabel(cb1,'Ice speed (m/yr)')

xlabel('x (m)');
ylabel('y (m)');
% tidal amplitude plot
ax2 = axes;

pcolor(ax2, x, y, KO_amp)
shading(ax2,'interp')
set(findobj(ax2,'Type','Surface'), ...
    'FaceAlpha',0.5, ...
    'EdgeColor','none')

colormap(ax2, turbo)
cb2 = colorbar(ax2,'eastoutside');
ylabel(cb2,'K1 + O1 amplitude (m)')

% Match axes exactly
set(ax2, ...
    'Position', ax1.Position, ...
    'Color','none', ...
    'XLim', xl, ...
    'YLim', yl, ...
    'DataAspectRatio', ax1.DataAspectRatio, ...
    'PlotBoxAspectRatio', ax1.PlotBoxAspectRatio, ...
    'XTick', [], ...
    'YTick', [], ...
    'Box','off')

linkaxes([ax1 ax2],'xy')

% Put ice-speed axes visually underneath
uistack(ax1,'bottom')

title(ax1,'CATS2008 K1 + O1 tidal amplitude over ice surface speed')


%% Now expand the constituent set to include all 10 constituents

clear;clc;

addpath('/Users/jeremywang/Documents/MATLAB/CATS2008')
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008';

latrange = -90.0:0.2:-60;
lonrange = -180:0.2:180;

[lon,lat] = meshgrid(lonrange,latrange);

% Check constituent list
[~,~,~,conList] = tmd_extract_HC(Model,lat(1),lon(1),'z');
conList_clean = lower(strtrim(cellstr(conList)));

% Extract harmonic amplitudes
[amp,~,~,~] = tmd_extract_HC(Model,lat,lon,'z');

% sum them up
ampsize = size(amp);
amp_sum = zeros(1, ampsize(2), ampsize(3));
for i=1:length(conList)
    amp_sum = amp_sum + amp(i,:,:);
end
amp_sum = squeeze(amp_sum);

% Convert CATS lon/lat grid to polar stereographic x/y used by CATS
[x,y] = ll2ps(lat,lon);

% add displacement
dx = 0;
x = x + dx;
figure
% ice speed plot
ax1 = axes;
measuresps('speed', 'log')

hold on
measuresps('gl','k')

axis(ax1,'equal')
axis(ax1,'manual')

xl = xlim(ax1);
yl = ylim(ax1);

colormap(ax1, parula)
cb1 = colorbar(ax1,'westoutside');
cb1.Ticks = [0 1 2 3];
cb1.TickLabels = {'1','10','100','1000'};
ylabel(cb1,'Ice speed (m/yr)')

xlabel('x (m)');
ylabel('y (m)');

% tidal amplitude plot
ax2 = axes;

pcolor(ax2, x, y, amp_sum)
shading(ax2,'interp')
set(findobj(ax2,'Type','Surface'), ...
    'FaceAlpha',0.5, ...
    'EdgeColor','none')

colormap(ax2, turbo)
cb2 = colorbar(ax2,'eastoutside');
ylabel(cb2,'Total tidal amplitude (m)')

% Match axes exactly
set(ax2, ...
    'Position', ax1.Position, ...
    'Color','none', ...
    'XLim', xl, ...
    'YLim', yl, ...
    'DataAspectRatio', ax1.DataAspectRatio, ...
    'PlotBoxAspectRatio', ax1.PlotBoxAspectRatio, ...
    'XTick', [], ...
    'YTick', [], ...
    'Box','off')

linkaxes([ax1 ax2],'xy')

% Put ice-speed axes visually underneath
uistack(ax1,'bottom')

title(ax1,'CATS2008 tidal amplitude over ice surface speed')

%% Calculate the beat frequency
clear;clc;
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008') 
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008'; 


latlim = [-76.0 -74.5];
lonlim = [-108 -103];

lats = latlim(1):0.5:latlim(2);
lons = lonlim(1):0.5:lonlim(2);

[~,~,~,conList] = tmd_extract_HC(Model,lats(1),lons(1),'z');
conList_clean = lower(strtrim(cellstr(conList)));

% start by calculating the beat frequency of just K1 and O1
% extract the periods
period = struct();
freq = struct();
for i = 1:length(conList_clean)
    con = conList_clean{i};
    % extract the period in hours
    [~,~,~,omega,~,~] = constit(con);
    T = 2*pi/omega * (1/(3600)); % converted to days
    
    freq.(con) = 1/T; % cycles per day 
    period.(con) = T; % days
end

% calculating beat frequency for each combination of constituents
beatFreq = struct();
beatPeriod = struct();
for i = 1:length(conList_clean)-1
    for j = i+1:length(conList_clean)

        con1 = conList_clean{i};
        con2 = conList_clean{j};

        name = sprintf('%s_%s', con1, con2);

        beatFreq.(name) = abs(freq.(con1) - freq.(con2)); % cycles/day
        beatPeriod.(name) = 1/beatFreq.(name);
    end
end

%% Extracting amplitudes of each constituent
extract_mean_amp('Thwaites Glacier');
beatPeriod.k1_o1
% around Thwaites the dominant tidal forcings are K1 and O1, meaning that
% we should use the beat period of 13.66 days

%% Generating map of daily tidal height for K1 and O1 in Antarctica
clear; clc;
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008')
addpath('/Users/jeremywang/Documents/MATLAB/measures_v3.1.2/measures')
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008';

latlim = [-78.0 -73];
lonlim = [-111 -99];

latrange = latlim(1):0.02:latlim(2);
lonrange = lonlim(1):0.02:lonlim(2);

[lon,lat] = meshgrid(lonrange,latrange);

% Check constituent list
conlist = extract_conlist(Model);

idx_k1 = find(strcmp(conlist,'k1'));
idx_o1 = find(strcmp(conlist,'o1'));

% time
t0 = datenum(2026,1,1,0,0,0);
time = t0 + (0:13);

% time = t0 + linspace(0,14,14);  % 14 snapshots over 14 days

% Extract 
height_total = NaN(length(time), length(latrange), length(lonrange));
for i=1:length(time)
    t = time(i);
    [h, ~] = tmd_tide_pred(Model, t, lat, lon, 'z', [idx_k1, idx_o1]);
    height_total(i,:,:) = h;
end

% Convert CATS lon/lat grid to polar stereographic x/y used by CATS
[x,y] = ll2ps(lat,lon);

for i=1:length(time)
    figure
    % ice speed plot
    ax1 = axes;
    mapzoomps('Thwaites Glacier')
    measuresps('speed','log')
    measuresps('gl', 'k')
    hold on
    
    axis(ax1,'equal')
    axis(ax1,'manual')
    
    xl = xlim(ax1);
    yl = ylim(ax1);
    
    colormap(ax1, parula)
    cb1 = colorbar(ax1,'westoutside');
    cb1.Ticks = [0 1 2 3];
    cb1.TickLabels = {'1','10','100','1000'};
    ylabel(cb1,'Ice speed (m/yr)')
    
    xlabel('x (m)');
    ylabel('y (m)');
    % tidal amplitude plot
    ax2 = axes;
    
    hmap = squeeze(height_total(i,:,:));
    pcolor(ax2, x, y, hmap)
    shading(ax2,'interp')
    set(findobj(ax2,'Type','Surface'), ...
        'FaceAlpha',0.5, ...
        'EdgeColor','none')
    
    colormap(ax2, turbo)
    cb2 = colorbar(ax2,'eastoutside');
    ylabel(cb2,'K1 + O1 amplitude (m)')
    maxabs = max(abs(height_total(:)), [], 'omitnan');
    clim = [-maxabs maxabs];

    
    % Match axes exactly
    set(ax2, ...
        'Position', ax1.Position, ...
        'Color','none', ...
        'XLim', xl, ...
        'YLim', yl, ...
        'DataAspectRatio', ax1.DataAspectRatio, ...
        'PlotBoxAspectRatio', ax1.PlotBoxAspectRatio, ...
        'XTick', [], ...
        'YTick', [], ...
        'Box','off')
    
    linkaxes([ax1 ax2],'xy')
    caxis(ax2, clim)

    % Put ice-speed axes visually underneath
    uistack(ax1,'bottom')
    
    title(ax1, sprintf('CATS2008 K1 + O1 tidal magnitude Day %.1f', time(i) - t0));
    figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/time_series_thwaites';
    exportgraphics(gcf, fullfile(figure_dir,sprintf('Thwaites_Day_%d.jpg', time(i) - t0)), ...
    'Resolution',300)
end
