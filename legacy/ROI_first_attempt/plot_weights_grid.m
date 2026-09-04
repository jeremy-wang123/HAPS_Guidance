function [] = plot_weights_grid(search_history, iteration)
%PLOT_WEIGHT_GRID creates two visualization plots of weighted grids
%   input struct that takes in the saved values and which iteration you
%   want

X_grid = search_history(iteration).X_grid;
Y_grid = search_history(iteration).Y_grid;
grounded_mask = search_history(iteration).grounded_mask;
weight_map = search_history(iteration).weight_map;

figure;

% Add the grounding line
measuresps('gl', 'k');
hold on;
scatter( ...
    X_grid(grounded_mask), ...
    Y_grid(grounded_mask), ...
    30, ...                           % marker size
    weight_map(grounded_mask), ...    % color values
    'filled', ...
    'MarkerEdgeColor', 'k' ...
);



mapzoomps('Thwaites Glacier');

cb = colorbar;
cb.Label.String = 'Weight score';

xlabel('Polar stereographic x (m)');
ylabel('Polar stereographic y (m)');
title(sprintf('Weight score over grounded ice in the Thwaites region at %d resolution', search_history(iteration).resolution));

fig_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/iceflowline_intersections';
exportgraphics(gcf, fullfile(fig_dir,sprintf('weight_grid_%d.jpg', iteration)), 'Resolution',300)
%%
x_grid = X_grid(1, :);   % First row
y_grid = Y_grid(:, 1).';   % First column

dx = x_grid(2) - x_grid(1);
dy = y_grid(2) - y_grid(1);

% Cell-edge coordinates, extending half a grid spacing around each center
x_edges = [x_grid - dx/2, x_grid(end) + dx/2];
y_edges = [y_grid - dy/2, y_grid(end) + dy/2];

[X_edges, Y_edges] = meshgrid(x_edges, y_edges);

% Pad the color matrix so pcolor does not discard the final data row/column
weight_padded = NaN(size(weight_map,1) + 1, ...
                    size(weight_map,2) + 1);

weight_padded(1:end-1, 1:end-1) = weight_map;

% Duplicate the boundary values into the padded row and column
weight_padded(end, 1:end-1) = weight_map(end,:);
weight_padded(1:end-1, end) = weight_map(:,end);
weight_padded(end,end) = weight_map(end,end);

figure;

measuresps('gl','k');

hold on; 
pcolor(X_edges, Y_edges, weight_padded);
shading flat;
axis image;
set(gca, 'YDir', 'normal');

% Plot the actual sample locations at the cell centers
valid_points = grounded_mask & isfinite(weight_map);

scatter( ...
    X_grid(valid_points), ...
    Y_grid(valid_points), ...
    20, ...
    'k', ...
    'filled');

mapzoomps('Thwaites Glacier');

cb = colorbar;
cb.Label.String = 'Weight score';

xlabel('Polar stereographic x (m)');
ylabel('Polar stereographic y (m)');
title(sprintf('Weight score over grounded ice in the Thwaites region at %d resolution', search_history(iteration).resolution));
exportgraphics(gcf, fullfile(fig_dir,sprintf('weight_dots_%d.jpg', iteration)), 'Resolution',300)
