%% Plotting ice stream boundaries

% Definition: A region in a grounded ice sheet in which the ice flows much faster than in the regions on either side.
% This means the ice stream must be behind the grounding line
% we can try a variety of velocity gradients to constrain flow that is
% "much faster"
% This lends itself to two variables: the distance scale considered, and
% the magnitude of the variability

%% Finding boundaries for ice streams
% Method: smooth the map from pixels, and then take the gradient at every
% point, take a minimum gradient


%% plot the velocity gradient 
clear; clc;

x = ncread('antarctica_ice_velocity_450m_v2.nc','x');
y = ncread(['antarctica_ice_velocity_450m_v' ...
    '2.nc'],'y');
u = ncread('antarctica_ice_velocity_450m_v2.nc','VX');
v = ncread('antarctica_ice_velocity_450m_v2.nc','VY');

u = u';
v = v';

% magnitude of the velocity vector
speed = hypot(u, v); 

% Grid spacing
dx = mean(diff(x));
dy = mean(diff(y));

% Gradient of speed
[dSdy, dSdx] = gradient(speed, dy, dx);

% Magnitude of velocity gradient
grad_speed = hypot(dSdx, dSdy);

%% histogram of distribution
figure
vals = grad_speed(isfinite(grad_speed));
histogram(log10(vals(vals>0)),100)
xlabel('log10 velocity gradient magnitude')
ylabel('Count')
title('Log distribution of velocity gradient')

%% plotting the velocity gradient
vals = grad_speed(isfinite(grad_speed));
percents = [90 95 97.5 99 99.5 99.9];
prctile(vals,percents);
thresholds = prctile(vals,percents);

fig_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/velocity_grad_plots';


for i=1:length(thresholds)
    figure;
    mapzoomps('Thwaites Glacier')
    xl = xlim;
    yl = ylim;
    
    
    measuresps('speed');
    hold on;
    colormap(parula);
    cb = colorbar;
    set(gca,'ColorScale');
    clim([1 1000]);
    cb.Label.String = 'Speed (m/yr)';
    measuresps('gl','k');

    contour(x,y,grad_speed,[thresholds(i) thresholds(i)],...
        'r','LineWidth',0.5);
    
    title(sprintf('Magnitude of velocity gradient for threshold of %.1f percent', percents(i)))
    xlim(xl);
    ylim(yl);
    xlabel('x (m)');
    ylabel('y (m)');
    % exportgraphics(gcf, fullfile(fig_dir,sprintf('velocity_grad_%d.jpg', percents(i)*10)), ...
    % 'Resolution',300)
end

%% Plotting strain rates

% Velocity gradients
[dudy, dudx] = gradient(u, dy, dx);
[dvdy, dvdx] = gradient(v, dy, dx);

% Strain-rate tensor components
exx = dudx;
eyy = dvdy;
exy = 0.5*(dudy + dvdx);

eff_strain = sqrt(exx.^2 + eyy.^2 + 2*exy.^2); % effective strain
shear_mag = abs(exy); % shear magnitude

plot_field = eff_strain;

vals = plot_field(isfinite(plot_field));
percentiles = [90 95 97.5 99 99.5 99.9];
thresholds = prctile(vals,percentiles);

for i=1:length(percentiles)
    fig = figure;
    mapzoomps('Thwaites Glacier')
    xl = xlim;
    yl = ylim;
    
    measuresps('speed')
    hold on;
    colormap(parula)
    colorbar
    set(gca,'ColorScale','log')
    clim([1 1000])

    measuresps('gl','k')

    contour(x,y,plot_field,[thresholds(i) thresholds(i)],'r','LineWidth',0.5)
    title(sprintf('Magnitude of effective strain for threshold of %.1f percent', percentiles(i)))
    
    xlim(xl)
    ylim(yl)
    xlabel('x (m)');
    ylabel('y (m)');
    
    fig_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/eff_strain_plots';

    exportgraphics(fig, fullfile(fig_dir,sprintf('eff_strain_%d.jpg', percentiles(i)*10)), ...
    'Resolution',300)
end
%% plotting grounded ice 

[xmesh, ymesh] = meshgrid(x, y);

grounded = isgrounded(xmesh, ymesh);

speed_plot = speed;

speed_grounded = speed_plot;
speed_grounded(~grounded) = NaN;

% Downsample factor
skip = 5;

x_ds = x(1:skip:end);
y_ds = y(1:skip:end);
speed_ds = speed_grounded(1:skip:end, 1:skip:end);

% recalculate gradient using grounded speed
[dSdy, dSdx] = gradient(speed_grounded, dy, dx);
grad_speed_grounded = hypot(dSdx, dSdy);

[xmesh_ds, ymesh_ds] = meshgrid(x_ds, y_ds);
percents = [90 95 97.5 99 99.5 99.9];

for i=1:length(percents)
    figure;
    pcolor(xmesh_ds, ymesh_ds, speed_ds);
    shading flat;
    axis image;
    set(gca,'YDir','normal');
    hold on
    
    colormap(parula);
    cb = colorbar;
    set(gca,'ColorScale');
    clim([1 1000]);
    cb.Label.String = 'Speed (m/yr)';
    
    mapzoomps('Thwaites Glacier')
    xl = xlim;
    yl = ylim;
    measuresps('gl','k')
    
    % velociy gradient lines
    vals = grad_speed_grounded(isfinite(grad_speed_grounded));
    
    prctile(vals,percents);
    thresholds = prctile(vals,percents);
    
    contour(x,y,grad_speed_grounded,[thresholds(i) thresholds(i)],'r','LineWidth',0.5)
    
    xlim(xl)
    ylim(yl)
    xlabel('x (m)');
    ylabel('y (m)');
    
    title(sprintf('Grounded ice velocity gradient for threshold of %.1f percent', percents(i)))
    
    fig_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/velocity_grad_plots';

    exportgraphics(gcf, fullfile(fig_dir,sprintf('grounded_velocity_grad_%d.jpg', percents(i)*10)), ...
    'Resolution',300)
end


%% Implementing a minimum pixel length

[xmesh, ymesh] = meshgrid(x, y);

grounded = isgrounded(xmesh, ymesh);

speed_plot = speed;

speed_grounded = speed_plot;
speed_grounded(~grounded) = NaN;

% Downsample factor
skip = 5;

x_ds = x(1:skip:end);
y_ds = y(1:skip:end);
speed_ds = speed_grounded(1:skip:end, 1:skip:end);

% recalculate gradient using grounded speed
[dSdy, dSdx] = gradient(speed_grounded, dy, dx);
grad_speed_grounded = hypot(dSdx, dSdy);

[xmesh_ds, ymesh_ds] = meshgrid(x_ds, y_ds);

% velociy gradient lines
vals = grad_speed_grounded(isfinite(grad_speed_grounded));
percents = [90 95 97.5 99 99.5 99.9];
prctile(vals,percents);
thresholds = prctile(vals,percents);


for i=1:length(percents)
    figure;
    pcolor(xmesh_ds, ymesh_ds, speed_ds);
    shading flat;
    axis image;
    set(gca,'YDir','normal');
    hold on
    
    colormap(parula);
    cb = colorbar;
    set(gca,'ColorScale');
    clim([1 1000]);
    cb.Label.String = 'Speed (m/yr)';
    
    mapzoomps('Thwaites Glacier')
    xl = xlim;
    yl = ylim;
    measuresps('gl','k')
    
    % Binary mask of high-gradient regions
    mask = grad_speed_grounded > thresholds(i);
    
    % Remove small connected components
    min_pix = 150;      % Adjust this value as needed
    mask = bwareaopen(mask, min_pix);
    
    % Plot contour of cleaned mask
    contour(x, y, double(mask), [1 1], ...
        'r', 'LineWidth', 0.5);
    
    xlim(xl)
    ylim(yl)
    xlabel('x (m)');
    ylabel('y (m)');
    
    title(sprintf('Grounded ice velocity gradient for threshold of %.1f percent', percents(i)))

    fig_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/velocity_grad_plots';
    
    % exportgraphics(gcf, fullfile(fig_dir,sprintf('cleaned_grounded_velocity_grad_%d.jpg', percents(i)*10)), ...
    % 'Resolution',300)
end





%% Localized gradient approach
clear; clc;

x = ncread('antarctica_ice_velocity_450m_v2.nc','x');
y = ncread(['antarctica_ice_velocity_450m_v' ...
    '2.nc'],'y');
u = ncread('antarctica_ice_velocity_450m_v2.nc','VX');
v = ncread('antarctica_ice_velocity_450m_v2.nc','VY');

u = u';
v = v';

% magnitude of the velocity vector
speed = hypot(u, v); 

% Grid spacing
dx = mean(diff(x));
dy = mean(diff(y));

[xmesh, ymesh] = meshgrid(x, y);

grounded = isgrounded(xmesh, ymesh);

speed_plot = speed;

speed_grounded = speed_plot;
speed_grounded(~grounded) = NaN;

% extract coordinate region of interest
fig = figure;
mapzoomps('Thwaites Glacier')
xl = xlim;
yl = ylim;
close(fig)

% Get x and y indices within the displayed region
ix = x >= xl(1) & x <= xl(2);
iy = y >= yl(1) & y <= yl(2);

% Subset coordinates
x_roi = x(ix);
y_roi = y(iy);

% Subset data
speed_roi = speed_grounded(iy, ix);

% Downsample factor
skip = 5;

x_ds = x_roi(1:skip:end);
y_ds = y_roi(1:skip:end);
speed_ds = speed_roi(1:skip:end, 1:skip:end);

% recalculate gradient using grounded speed
[dSdy, dSdx] = gradient(speed_roi, dy, dx);
grad_speed_grounded = hypot(dSdx, dSdy);

[xmesh_ds, ymesh_ds] = meshgrid(x_ds, y_ds);

% establish percentile distributions
percents = [90 95 97.5 99 99.5 99.9];

% establish a window in which each pixel is compared to
window = [100 100];

% local percentile rank of each pixel
local_rank = nlfilter(grad_speed_grounded, window, @(z) local_percentile_rank(z));

% fig_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/velocity_grad_plots';

mask = local_rank > percents(3)/100;

% remove tiny noisy patches
min_pix = 0; % determines the minimum pixels to include a line
mask = bwareaopen(mask, min_pix);

figure;
% first plotting the downsized velocity profile
pcolor(xmesh_ds, ymesh_ds, speed_ds);
shading flat;
axis image
set(gca,'YDir','normal')
hold on

colormap(parula);
cb = colorbar;
set(gca,'ColorScale');
clim([1 1000]);
cb.Label.String = 'Speed (m/yr)';

% selecting for the Thwaites region
mapzoomps('Thwaites Glacier')
xl = xlim;
yl = ylim;
measuresps('gl','k'); % grounding line

% plotting velocity gradient contours
contour(x_roi,y_roi,double(mask),[1 1], ...
    'r','LineWidth',0.5);

title(sprintf('Local velocity-gradient percentile threshold of %.1f percent', percents(2)))
xlim(xl);
ylim(yl);
xlabel('x (m)');
ylabel('y (m)');

% this function takes in an individual pixel and outputs the percentile
% distribution the pixel is, compared to its surrounding neighbors
function r = local_percentile_rank(z)
    center = z(ceil(numel(z)/2));
    vals = z(isfinite(z));

    if ~isfinite(center) || isempty(vals)
        r = NaN;
    else
        r = mean(vals <= center);
    end
end

%% Using antbounds and ice flowlines package

%%% Reading data
x = ncread('antarctica_ice_velocity_450m_v2.nc','x');
y = ncread(['antarctica_ice_velocity_450m_v' ...
    '2.nc'],'y');
u = ncread('antarctica_ice_velocity_450m_v2.nc','VX');
v = ncread('antarctica_ice_velocity_450m_v2.nc','VY');

u = u';
v = v';

% magnitude of the velocity vector
speed = hypot(u, v); 

[xmesh, ymesh] = meshgrid(x, y);

%%% Eliminating data points not on grounded ice
grounded = isgrounded(xmesh, ymesh);
speed_plot = speed;
speed_grounded = speed_plot;
speed_grounded(~grounded) = NaN;

%%% extract coordinate region of interest
fig = figure;
mapzoomps('Thwaites Glacier')
xl = xlim;
yl = ylim;
close(fig)

% Get x and y indices within the displayed region
ix = x >= xl(1) & x <= xl(2);
iy = y >= yl(1) & y <= yl(2);

% Subset coordinates
x_roi = x(ix);
y_roi = y(iy);

% Subset data
speed_roi = speed_grounded(iy, ix);

figure;
pcolor(x_roi, y_roi, speed_roi);
shading flat;
axis image;
set(gca,'YDir','normal');

colormap(parula);
cb = colorbar;
set(gca,'ColorScale');
clim([1 1000]);
cb.Label.String = 'Speed (m/yr)';

hold on;

% plotting ice flowlines

[xstart,ystart] = psgrid('thwaites glacier',500,30, 'xy');

flowline(xstart,ystart,'plotxy', 'color', 'green')

% plotting the ice shelf
[xgrid,ygrid] = psgrid('Thwaites Glacier',700,10,'xy');
axis image
antbounds('gl','k')
antbounds('coast','k')
antbounds('shelves','k')
scalebarps

shelf = isiceshelf(xgrid,ygrid);
plot(xgrid(shelf),ygrid(shelf),'kx')

[wx,wy] = antbounds_data('Thwaites','xy');
plot(wx,wy,'r','linewidth',2)

thwaites = inpolygon(xgrid,ygrid,wx,wy);
plot(xgrid(thwaites),ygrid(thwaites),'ro')

% grounding line
mapzoomps('Thwaites Glacier')
xl = xlim;
yl = ylim;
measuresps('gl','k');

% labels
xlabel('x (m)');
ylabel('y (m)');
xlim(xl);
ylim(yl);
