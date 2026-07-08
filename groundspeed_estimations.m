%% Plotting estimations for ground speed
t = 0:0.1:20; % days
speed = 20; % m/s (estimation)
v = speed*(60*60*24)/(1000); % km/day 

figure;
plot(t, v*t);
xlabel('Time (days)');
ylabel('Distance (km)');
title('Flight Range');