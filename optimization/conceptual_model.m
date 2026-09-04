%% Conceptual optimization figure
% Tile-time decision space over one day

rng(4);   % Reproducible random pattern

n_tiles = 16;
total_days = 1;

figure;
hold on;


%% Grid boundaries

% Vertical tile boundaries
for i = 0:n_tiles
    xline(i + 0.5, ...
        'Color',[0.8 0.8 0.8], ...
        'LineWidth',0.75, ...
        'HandleVisibility','off');
end

% Horizontal time guides every 3 hours
time_ticks_hr = 0:3:24;
time_ticks_day = time_ticks_hr/24;

for t = time_ticks_day
    yline(t, ...
        'Color',[0.85 0.85 0.85], ...
        'LineWidth',0.75, ...
        'HandleVisibility','off');
end


%% Generate somewhat-random visits

for tile = 1:n_tiles

    % First visit occurs sometime in first 8 hours
    first_visit = 8*rand;       % [hr]

    % Second visit roughly 12 hours later,
    % with +/- 2 hours of random variation
    second_visit = first_visit + 12 + 4*(rand - 0.5);

    visit_times_hr = first_visit;

    % Only include second visit if it occurs within the day
    if second_visit <= 24
        visit_times_hr(end+1) = second_visit;
    end


    %% Plot visits

    for t_hr = visit_times_hr

        t_day = t_hr/24;

        % Line spans width of tile column
        plot( ...
            [tile-0.45 tile+0.45], ...
            [t_day t_day], ...
            'k-', ...
            'LineWidth',1.5, ...
            'HandleVisibility','off');

    end

end


%% Axes

xlim([0.5 n_tiles+0.5]);
ylim([0 1]);

set(gca, ...
    'XTick',1:n_tiles, ...
    'YTick',time_ticks_day, ...
    'YTickLabel',string(time_ticks_hr), ...
    'XAxisLocation','top', ...
    'YDir','reverse');

xlabel('Tile Number');
ylabel('Time (hours)');

title('Tile-Time Decision Space');

box on;

set(gca,'FontSize',11);


figure_dir = fullfile('/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/optimization_figures');

% Create folder if it does not exist
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

exportgraphics(gcf, ...
   fullfile(figure_dir, 'tile_time_decisionspace.png'), ...
   'Resolution', 300);