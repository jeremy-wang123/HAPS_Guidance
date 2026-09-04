%% Tidal Constituents Animation
%% Time
t = linspace(0,72,1500);   % hours

%% Tidal periods [hours]
% extract the periods
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008'; 
conList_clean = extract_conlist(Model);
period = struct();
freq = struct();
for i = 1:length(conList_clean)
    con = conList_clean{i};
    % extract the period in hours
    [~,~,~,omega,~,~] = constit(con);
    T = 2*pi/omega * (1/(3600)); 
    
    freq.(con) = 1/T; % cycles per hour 
    period.(con) = T; % hour
end

%% Amplitudes [m]
mean_amps = extract_mean_amp('Thwaites Glacier');
%% Signals

M2 = mean_amps.m2 .* sin(2*pi*t ./ period.m2);
S2 = mean_amps.s2 .* sin(2*pi*t ./ period.s2);
K1 = mean_amps.k1 .* sin(2*pi*t ./ period.k1);
O1 = mean_amps.o1 .* sin(2*pi*t ./ period.o1);
K2 = mean_amps.k2 .* sin(2*pi*t ./ period.k2);

composite = M2 + S2 + K1 + O1 + K2;

%% colors
dark_green = [0.00 0.45 0.25];
purple     = [0.50 0.15 0.70];
dark_blue  = [0.05 0.25 0.65];
orange     = [0.95 0.40 0.05];
yellow     = [0.90 0.70 0.00];

color.M2 = dark_green;
color.S2 = purple;
color.K1 = dark_blue;
color.O1 = orange;
color.K2 = yellow;
%% Create figure
fig = figure( ...
    'Position',[100 100 1100 750], ...
    'Color','w');

tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% ---------------- TOP PLOT ----------------
ax1 = nexttile;

plot(t,M2,'Color',color.M2, ...
    'LineWidth',1.5, ...
    'DisplayName','M_2 (12.42 h)');
hold on

plot(t,S2,'Color', color.S2, ...
    'LineWidth',1.5, ...
    'DisplayName','S_2 (12.00 h)');

plot(t,K1,'Color', color.K1, ...
    'LineWidth',1.5, ...
    'DisplayName','K_1 (23.93 h)');

plot(t,O1,'Color', color.O1,...
    'LineWidth',1.5, ...
    'DisplayName','O_1 (25.82 h)');

plot(t,K2,'Color', color.K2, ...
    'LineWidth',1.5, ...
    'DisplayName','K_2 (11.97 h)');


xlim([0 72])
ylim([-0.5 0.5])

xlabel('Time (hours)')
ylabel('Height (m)')

title('Individual Tidal Constituents', ...
    'FontSize',16, ...
    'FontWeight','bold')

legend('Location','northeast')

grid on
box on

% ---------------- BOTTOM PLOT ----------------
ax2 = nexttile;

plot(t,composite,'k','LineWidth',2);
hold on

yline(0,'k--');

xlim([0 72])
ylim([-1 1])

xlabel('Time (hours)')
ylabel('Height (m)')

title('Composite Tide (Sum of All Constituents)', ...
    'FontSize',16, ...
    'FontWeight','bold')

grid on
box on


%% Dominant Constituent Figure
fig = figure( ...
    'Position',[100 100 1100 750], ...
    'Color','w');

tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% ---------------- TOP PLOT ----------------
ax1 = nexttile;

plot(t,M2,'Color',color.M2, ...
    'LineWidth',0.5, ...
    'DisplayName','M_2 (12.42 h)');
hold on

plot(t,S2,'Color', color.S2, ...
    'LineWidth',0.5, ...
    'DisplayName','S_2 (12.00 h)');

plot(t,K1,'Color', color.K1, ...
    'LineWidth',3, ...
    'DisplayName','K_1 (23.93 h)');

plot(t,O1,'Color', color.O1,...
    'LineWidth',0.5, ...
    'DisplayName','O_1 (25.82 h)');

plot(t,K2,'Color', color.K2, ...
    'LineWidth',0.5, ...
    'DisplayName','K_2 (11.97 h)');

xline(period.k1/2, 'r--', 'LineWidth', 3)

xlim([0 72])
ylim([-0.5 0.5])

xlabel('Time (hours)')
ylabel('Height (m)')

title('Individual Tidal Constituents', ...
    'FontSize',16, ...
    'FontWeight','bold')

legend('Location','northeast')

grid on
box on

% ---------------- BOTTOM PLOT ----------------
ax2 = nexttile;

plot(t,composite,'k','LineWidth',2);
hold on

yline(0,'k--');

xlim([0 72])
ylim([-1 1])

xlabel('Time (hours)')
ylabel('Height (m)')

title('Composite Tide (Sum of All Constituents)', ...
    'FontSize',16, ...
    'FontWeight','bold')

grid on
box on

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/animations';

exportgraphics(gcf, fullfile(figure_dir,'dominant_const.jpg'), ...
'Resolution',300);

%% K1 and O1 plot combined
t = linspace(0, 25*24, 1000);
t_days = t/24; 
K1 = mean_amps.k1 .* sin(2*pi*t ./ period.k1);
O1 = mean_amps.o1 .* sin(2*pi*t ./ period.o1);

K1_O1 = K1 + O1;

figure('Position',[100 100 1000 800]);

% ---------------- K1 ----------------
subplot(3,1,1)

plot(t_days, K1, ...
    'Color',[0.05 0.25 0.65], ...
    'LineWidth',2);

yline(0,'k--');

xlabel('Time (Days)')
ylabel('Height (m)')
title(sprintf('K_1 Constituent (T = %.2f h)', period.k1))

grid on
box on


% ---------------- O1 ----------------
subplot(3,1,2)

plot(t_days, O1, ...
    'Color',[0.95 0.40 0.05], ...
    'LineWidth',2);

yline(0,'k--');
xlabel('Time (Days)')
ylabel('Height (m)')
title(sprintf('O_1 Constituent (T = %.2f h)', period.o1))

grid on
box on


% ---------------- K1 + O1 ----------------
subplot(3,1,3)

% Beat period
f_K1 = freq.k1;
f_O1 = freq.o1;

T_beat_hours = 1 / abs(f_K1 - f_O1);
T_beat_days  = T_beat_hours / 24;

% Offset so first line is at destructive interference
offset = T_beat_days / 2;

beat_times = offset:T_beat_days:max(t_days);


plot(t_days, K1_O1, ...
    'Color',[0.50 0.15 0.70], ...
    'LineWidth',2);

for tb = beat_times
    xline(tb, '--r', ...
            'LineWidth',3, ...
            'HandleVisibility','off');
end
yline(0,'k--');


xlabel('Time (Days)')
ylabel('Height (m)')
title('K1 + O1')
grid on
box on


figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/animations';

exportgraphics(gcf, fullfile(figure_dir,'beat_freq.jpg'), ...
'Resolution',300);

%% Individual Figures
% Beat period
f_K1 = freq.k1;
f_O1 = freq.o1;

T_beat_hours = 1 / abs(f_K1 - f_O1);
T_beat_days  = T_beat_hours / 24;

% Offset so first line is at destructive interference
offset = T_beat_days / 2;

beat_times = offset:T_beat_days:max(t_days);


plot(t_days, K1_O1, ...
    'Color',[0.50 0.15 0.70], ...
    'LineWidth',2);

for tb = beat_times
    xline(tb, '--r', ...
            'LineWidth',3, ...
            'HandleVisibility','off');
end
yline(0,'k--');


xlabel('Time (Days)')
ylabel('Height (m)')
title('Minimum Campaign Period')
grid on
box on
figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/animations';

exportgraphics(gcf, fullfile(figure_dir,'individual_beat_freq.jpg'), ...
'Resolution',300);

%%
fig = figure( ...
    'Position',[100 100 1100 750], ...
    'Color','w');

tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


plot(t,M2,'Color',color.M2, ...
    'LineWidth',0.5, ...
    'DisplayName','M_2 (12.42 h)');
hold on

plot(t,S2,'Color', color.S2, ...
    'LineWidth',0.5, ...
    'DisplayName','S_2 (12.00 h)');

plot(t,K1,'Color', color.K1, ...
    'LineWidth',3, ...
    'DisplayName','K_1 (23.93 h)');

plot(t,O1,'Color', color.O1,...
    'LineWidth',0.5, ...
    'DisplayName','O_1 (25.82 h)');

plot(t,K2,'Color', color.K2, ...
    'LineWidth',0.5, ...
    'DisplayName','K_2 (11.97 h)');

xline(period.k1/2, 'r--', 'LineWidth', 3,'HandleVisibility','off')
xline(period.k1, 'r--', 'LineWidth', 3, 'HandleVisibility','off')

xlim([0 72])
ylim([-0.5 0.5])

xlabel('Time (hours)')
ylabel('Height (m)')

title('Maximum Sampling Period', ...
    'FontSize',16, ...
    'FontWeight','bold')

legend('Location','northeast')

grid on
box on

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/animations';

exportgraphics(gcf, fullfile(figure_dir,'individual_dom_const.jpg'), ...
'Resolution',300);

%% K1 and O1 animated plot

t = linspace(0, 25*24, 1000);   % hours
t_days = t/24;

K1 = mean_amps.k1 .* sin(2*pi*t ./ period.k1);
O1 = mean_amps.o1 .* sin(2*pi*t ./ period.o1);

K1_O1 = K1 + O1;

% Figure
fig = figure( ...
    'Position',[100 100 1000 800], ...
    'Color','w');
% ---------------- K1 ----------------
ax1 = subplot(3,1,1);

h1 = plot(t_days(1), K1(1), ...
    'Color',[0.05 0.25 0.65], ...
    'LineWidth',2);

hold on
yline(0,'k--');

xlim([t_days(1) t_days(end)])
ylim([-1.1*max(abs(K1)) 1.1*max(abs(K1))])

xlabel('Time (Days)')
ylabel('Height (m)')
title(sprintf('K_1 Constituent (T = %.2f h)', period.k1))

grid on
box on


% ---------------- O1 ----------------
ax2 = subplot(3,1,2);

h2 = plot(t_days(1), O1(1), ...
    'Color',[0.95 0.40 0.05], ...
    'LineWidth',2);

hold on
yline(0,'k--');

xlim([t_days(1) t_days(end)])
ylim([-1.1*max(abs(O1)) 1.1*max(abs(O1))])

xlabel('Time (Days)')
ylabel('Height (m)')
title(sprintf('O_1 Constituent (T = %.2f h)', period.o1))

grid on
box on


% ---------------- K1 + O1 ----------------
ax3 = subplot(3,1,3);

h3 = plot(t_days(1), K1_O1(1), ...
    'Color',[0.50 0.15 0.70], ...
    'LineWidth',2);

hold on
yline(0,'k--');

xlim([t_days(1) t_days(end)])
ylim([-1.1*max(abs(K1_O1)) 1.1*max(abs(K1_O1))])

xlabel('Time (Days)')
ylabel('Height (m)')
title('K_1 + O_1')

grid on
box on


% Moving dots
dot1 = plot(ax1, t_days(1), K1(1), 'o', ...
    'Color',[0.05 0.25 0.65], ...
    'MarkerFaceColor',[0.05 0.25 0.65], ...
    'MarkerSize',7);

dot2 = plot(ax2, t_days(1), O1(1), 'o', ...
    'Color',[0.95 0.40 0.05], ...
    'MarkerFaceColor',[0.95 0.40 0.05], ...
    'MarkerSize',7);

dot3 = plot(ax3, t_days(1), K1_O1(1), 'o', ...
    'Color',[0.50 0.15 0.70], ...
    'MarkerFaceColor',[0.50 0.15 0.70], ...
    'MarkerSize',7);


% GIF settings
figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/animations';

if ~exist(figure_dir,'dir')
    mkdir(figure_dir)
end

fname = fullfile(figure_dir,'K1_O1_animation');

framerate = 30;

% Use fewer animation frames than total data points
frame_idx = round(linspace(1,length(t),300));


% Animate
for i = 1:length(frame_idx)

    idx = frame_idx(i);

    % Draw the curves progressively
    h1.XData = t_days(1:idx);
    h1.YData = K1(1:idx);

    h2.XData = t_days(1:idx);
    h2.YData = O1(1:idx);

    h3.XData = t_days(1:idx);
    h3.YData = K1_O1(1:idx);

    % Move dots
    dot1.XData = t_days(idx);
    dot1.YData = K1(idx);

    dot2.XData = t_days(idx);
    dot2.YData = O1(idx);

    dot3.XData = t_days(idx);
    dot3.YData = K1_O1(idx);

    drawnow

    % Save frame
    if i == 1
        gifanim(fname, framerate, fig, 1);
    else
        gifanim(fname, framerate, fig, 0);
    end

end

%% -------------------------------------------
% Animation objects
%% -------------------------------------------

% Vertical moving bars
bar1 = xline(ax1,0, ...
    'k:', ...
    'LineWidth',2,  ...
    'HandleVisibility','off');

bar2 = xline(ax2,0, ...
    'k:', ...
    'LineWidth',2);

% Moving dot on composite tide
dot = plot(ax2, ...
    t(1),composite(1), ...
    'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',5);

% Moving dots on individual constituents

dot_M2 = plot(ax1, t(1), M2(1), 'o', ...
    'MarkerFaceColor',color.M2, ...
    'MarkerSize',6, ...
    'HandleVisibility','off');

dot_S2 = plot(ax1, t(1), S2(1), 'bo', ...
    'MarkerFaceColor',color.S2, ...
    'MarkerSize',6, ...
    'HandleVisibility','off');

dot_K1 = plot(ax1, t(1), K1(1), 'go', ...
    'MarkerFaceColor',color.K1, ...
    'MarkerSize',6, ...
    'HandleVisibility','off');

dot_O1 = plot(ax1, t(1), O1(1), 'o', ...
    'Color',[1 0.25 0], ...
    'MarkerFaceColor', color.O1, ...
    'MarkerSize',6, ...
    'HandleVisibility','off');

dot_K2 = plot(ax1, t(1), K2(1), 'o', ...
    'Color',[1 0.6 0], ...
    'MarkerFaceColor',color.K2, ...
    'MarkerSize',6, ...
    'HandleVisibility','off');

% Time label above upper plot
time_top = text(ax1, ...
    0,1.55, ...
    '0.0 h', ...
    'Color',[0.55 0 0.8], ...
    'FontWeight','bold', ...
    'FontSize',12, ...
    'HorizontalAlignment','center');

% Bottom time label
time_box = annotation(fig, ...
    'textbox', ...
    [0.40 0.005 0.20 0.04], ...
    'String','Time = 0.0 hours', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontWeight','bold', ...
    'FontSize',12, ...
    'BackgroundColor','w', ...
    'EdgeColor','k', ...
    'FitBoxToText','on');

%% -------------------------------------------
% GIF settings
%% -------------------------------------------

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/animations';
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

fname = fullfile(figure_dir, 'tidal_constituents');
framerate = 12;

nframes = 240;

frame_times = linspace(0,72,nframes);

%% -------------------------------------------
% Animate
%% -------------------------------------------
for i = 1:nframes

    current_time = frame_times(i);

    % Find closest point
    [~,idx] = min(abs(t-current_time));

    % Move vertical bars
    bar1.Value = current_time;
    bar2.Value = current_time;

    % Move dots on individual constituents
    dot_M2.XData = current_time;
    dot_M2.YData = M2(idx);

    dot_S2.XData = current_time;
    dot_S2.YData = S2(idx);

    dot_K1.XData = current_time;
    dot_K1.YData = K1(idx);

    dot_O1.XData = current_time;
    dot_O1.YData = O1(idx);

    dot_K2.XData = current_time;
    dot_K2.YData = K2(idx);

    % Move composite point
    dot.XData = current_time;
    dot.YData = composite(idx);

    % Update time labels
    time_top.Position(1) = current_time;
    time_top.String = sprintf('%.1f h',current_time);

    time_box.String = ...
        sprintf('Time = %.1f hours',current_time);

    drawnow

    % Write GIF
    if i == 1
        gifanim(fname,framerate,fig,1);
    else
        gifanim(fname,framerate,fig,0);
    end

end

%% viscous beam bending
% Number of animation frames
nFrames = 120;

% Go through one complete tidal cycle
phases_deg = linspace(0, 360, nFrames);

% Twilight colormap
cmap = slanCM('twilight_s',256);
nColors = size(cmap,1);

%% -------------------------------------------------
% Determine fixed y limits before animation
%% -------------------------------------------------

% Maximum possible displacement envelope
w_amp = abs(W);

w_max = max(w_amp);

vertical_padding = 0.15*w_max;

y_limits = [ ...
    -w_max - vertical_padding, ...
     w_max + vertical_padding];

%% -------------------------------------------------
% Create figure
%% -------------------------------------------------

fig = figure( ...
    'Position',[100 100 1000 650], ...
    'Color','w');

ax = axes(fig);
hold(ax,'on');

% Initial phase
phase_deg = phases_deg(1);
t_current = (phase_deg/360)*T;

w = imag(W .* exp(1i*omega*t_current));

% Initial color
color_index = 1;
line_color = cmap(color_index,:);

% Plot displacement
h = plot(x_km, w, ...
    'Color',line_color, ...
    'LineWidth',3);

% Grounding line / zero displacement
yline(0,'k--', ...
    'HandleVisibility','off');

%% Formatting
xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');

title('Thwaites Glacier Tidal Deflection');

xlim([0 L/1000]);
ylim(y_limits);

grid on;
box on;

set(gca, ...
    'FontSize',12, ...
    'LineWidth',1);

%% Colorbar
colormap(gca,cmap);
clim([0 360]);

cb = colorbar;
cb.Label.String = 'Tidal phase (degrees)';

cb.Ticks = 0:90:360;
cb.TickLabels = { ...
    '0^\circ', ...
    '90^\circ', ...
    '180^\circ', ...
    '270^\circ', ...
    '360^\circ'};

%% Moving phase label
phase_text = text( ...
    0.03,0.92, ...
    'Phase = 0^\circ', ...
    'Units','normalized', ...
    'FontSize',14, ...
    'FontWeight','bold');

%% -------------------------------------------------
% GIF settings
%% -------------------------------------------------

save_dir = '/Users/jeremywang/Desktop/GIFs';

if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

fname = fullfile(save_dir,'tidal_deflection');

framerate = 25;

%% -------------------------------------------------
% Animate
%% -------------------------------------------------

for j = 1:nFrames

    % Current tidal phase
    phase_deg = phases_deg(j);

    % Convert phase to time
    t_current = (phase_deg/360)*T;

    % Calculate instantaneous displacement
    w = imag(W .* exp(1i*omega*t_current));

    % Select color corresponding to phase
    color_position = phase_deg/360;

    color_index = 1 + ...
        round(color_position*(nColors-1));

    line_color = cmap(color_index,:);

    % Update displacement
    h.YData = w;

    % Update line color
    h.Color = line_color;

    % Update phase label
    phase_text.String = ...
        sprintf('Phase = %.0f^\\circ',phase_deg);

    drawnow limitrate

    %% Write GIF
    if j == 1
        gifanim(fname,framerate,fig,1);
    else
        gifanim(fname,framerate,fig,0);
    end

end
%%
function [] = gifanim(fname, framerate, fighandl, loopind)

% Capture frame
frame = getframe(fighandl);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

if loopind == 1

    imwrite(imind,cm,[fname '.gif'],'gif', ...
        'LoopCount',inf, ...
        'DelayTime',1/framerate);

else

    imwrite(imind,cm,[fname '.gif'],'gif', ...
        'WriteMode','append', ...
        'DelayTime',1/framerate);

end

end