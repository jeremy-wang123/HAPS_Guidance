function [] = plot_line_trajectory(distance, theta, x0, y0)
%PLOT TRAJECTORY plots line trajectory of candidate flight
% limited to the Thwaites region
x1 = x0 + distance*cosd(theta);
y1 = y0 + distance*sind(theta);

figure
mapzoomps('Thwaites Glacier')

% Background speed map
measuresps('speed')
hold on

% Grounding line
measuresps('gl','k')

% Flowlines
[xstart,ystart] = psgrid('thwaites glacier',500,30,'xy');
flowline(xstart,ystart,'plotxy','color','g')

% Candidate flight line
plot([x0 x1], [y0 y1], ...
    'r-', 'LineWidth', 4)

% mark endpoints
plot(x0,y0,'ro','MarkerFaceColor','r','MarkerSize',8)
plot(x1,y1,'ro','MarkerFaceColor','r','MarkerSize',8)

xlabel('x (m)')
ylabel('y (m)')
title('Candidate Straight Flight Line Over Thwaites')
