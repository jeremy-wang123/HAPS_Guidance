%% Preliminary plots

%% Quiver and color plot of ice velocities
figure
measuresps('speed','log','alpha',0.7)
colorbar
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
skip = 200;   

xq = x(1:skip:end);
yq = y(1:skip:end);
vxq = vx(1:skip:end,1:skip:end);
vyq = vy(1:skip:end,1:skip:end);

[X,Y] = meshgrid(xq,yq);

h = quiver(X,Y,vxq',vyq','k','LineWidth',0.8);
xlim(xl)
ylim(yl)
set(h,'AutoScaleFactor',1)
axis equal


