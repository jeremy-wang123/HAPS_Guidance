%% Overlay the maximum distance of flight path in one week

%%  calculate the flight distance for straight line flight
campaign = 7; % 7 day campaign
speed = 20; % m/s (estimation)
v = speed*(60*60*24); % m/day 

max_distance = v*campaign/2; % m
%% Simple plot of Thwaites region
measuresps('speed')
mapzoomps('Thwaites Glacier')
hold on;
measuresps('gl', 'k')

colormap(parula);
cb = colorbar;
set(gca,'ColorScale');
clim([1 1000]);
cb.Label.String = 'Speed (m/yr)';

xlabel('x (m)');
ylabel('y (m)');

x0 = -1.5e6;
y0 = -4e5;
hold on
plot([x0, x0+max_distance], [y0, y0], 'r-', 'LineWidth', 3)


%% Basic idea:
% we follow the grounding line of the ice sheet (probably starting near
% Thwaites) and use each point as a starting point for the drone takeoff.
% We then vary the direction of flight by a fraction of a degree,
% effectively simulating flight in different directions. Within a given
% trajectory, we determine how many ice flowline intersections there are,
% and we weight the ice flowlines according to relevancy (this can be
% determined by the proximinity to the grounding line). Then we can choose
% what options are the best

% This method will certainly be very computationally expensive, so not sure
% we should go ahead with it.

x0 = -1.5e6;
y0 = -6e5;
distance = 3e5;
theta = 90;

plot_line_trajectory(distance,theta,x0,y0)


%% Creating different potential flight paths in Thwaites
% each scheduled flight has two variables: starting location, and direction
% of flight, assuming a straight line

% this is obviously computationally infeasible, so we need ot figure out
% some ways of eliminating starting points and/or directions

% 1) We can elimate entirely any point on the map that is not on grounded
% ice. The aircraft will need to take off on a runway, which means anything
% on water will not work
% 2) We can also sweep only from 0 to 180 degrees, because the coverage
% would be the same. This halves our compute time, but also makes sure we
% do not get 2 local maximums (since the same line would be found in 2
% different locations)

%% Iterative coarse-to-fine weight search

distance = 3e5;

% Initial Thwaites search region
figure;
mapzoomps('Thwaites Glacier');
x_limits = xlim;
y_limits = ylim;
close;

% Search settings
grid_resolution = 100e3;       % initial spacing: 100 km
refinement_factor = 5;         % divide spacing by 5 each iteration
minimum_resolution = 1e3;      % stop at 1 km resolution
max_iterations = 10;
angle_spacing = 45; % 45 degrees between angles being tested

% Store results from each refinement level
search_history = struct([]);

for iteration = 1:max_iterations
    fprintf('\nIteration %d\n', iteration);
    fprintf('Grid resolution: %.2f km\n', grid_resolution/1e3);
    
    fprintf('x Limits: [%.1f, %.1f] m\n', x_limits(1), x_limits(2));
    fprintf('y Limits: [%.1f, %.1f] m\n', y_limits(1), y_limits(2));

    %%% Construct current grid

    x_grid = x_limits(1):grid_resolution:x_limits(2);
    y_grid = y_limits(1):grid_resolution:y_limits(2);

    [X_grid, Y_grid] = meshgrid(x_grid, y_grid);

    grounded_mask = isgrounded(X_grid, Y_grid);
    weight_map = NaN(size(X_grid));

    grounded_idx = find(grounded_mask);
    N = numel(grounded_idx);

    if N == 0
        error('No grounded grid points exist in the current search region.');
    end

    %%% Evaluate weights

    h = waitbar(0, sprintf('Iteration %d...', iteration));
    tStart = tic;

    for k = 1:N

        idx = grounded_idx(k);

        x_start = X_grid(idx);
        y_start = Y_grid(idx);

        weight_map(idx) = cumscore( ...
            distance, ...
            angle_spacing, ...
            x_start, ...
            y_start);

        if mod(k,10) == 0 || k == N

            elapsed = toc(tStart);
            remaining = elapsed/k*(N-k);

            waitbar(k/N, h, sprintf( ...
                ['Iteration %d\n' ...
                 '%.1f%% complete\n' ...
                 'Elapsed: %.1f min\n' ...
                 'Remaining: %.1f min'], ...
                iteration, ...
                100*k/N, ...
                elapsed/60, ...
                remaining/60));
        end
    end

    close(h);

    %%% Find current maximum

    [max_weight, linear_idx] = max( ...
        weight_map(:), ...
        [], ...
        'omitnan');

    x_best = X_grid(linear_idx);
    y_best = Y_grid(linear_idx);

    [best_row, best_col] = ind2sub( ...
        size(weight_map), ...
        linear_idx);

    fprintf('Best weight: %.6f\n', max_weight);
    fprintf('Best point: x = %.1f m, y = %.1f m\n', ...
        x_best, y_best);

    %%% Save this iteration

    search_history(iteration).resolution = grid_resolution;
    search_history(iteration).x_limits = x_limits;
    search_history(iteration).y_limits = y_limits;
    search_history(iteration).X_grid = X_grid;
    search_history(iteration).Y_grid = Y_grid;
    search_history(iteration).grounded_mask = grounded_mask;
    search_history(iteration).weight_map = weight_map;
    search_history(iteration).max_weight = max_weight;
    search_history(iteration).x_best = x_best;
    search_history(iteration).y_best = y_best;
    search_history(iteration).best_row = best_row;
    search_history(iteration).best_col = best_col;

    %%%  Check stopping condition

    if grid_resolution <= minimum_resolution
        fprintf('\nMinimum resolution reached.\n');
        break
    end

    %%%  Refine around the best point

    coarse_resolution = grid_resolution;
    grid_resolution = grid_resolution / refinement_factor;


    % Search one full coarse-grid spacing around the current best point
    search_half_width = coarse_resolution;

    x_limits = [
        x_best - search_half_width, ...
        x_best + search_half_width
    ];

    y_limits = [
        y_best - search_half_width, ...
        y_best + search_half_width
    ];
    
end

%% Final result

final_x = x_best;
final_y = y_best;
final_weight = max_weight;
final_resolution = grid_resolution;

fprintf('\nFinal selected point:\n');
fprintf('x = %.3f m\n', final_x);
fprintf('y = %.3f m\n', final_y);
fprintf('Weight = %.6f\n', final_weight);
fprintf('Spatial resolution = %.3f km\n', ...
    final_resolution/1e3);
 
% finer iterations through specific angles
angles = 0:5:175;
angle_weights = NaN(size(angles));

for i=1:length(angles)
    theta = angles(i);

    angle_weights(i) = extract_weights(distance,theta,final_x,final_y);
end

% Find best angle
[max_weight, best_idx] = max(angle_weights);

best_angle = angles(best_idx);

fprintf('\nOptimal location:\n');
fprintf('x = %.1f m\n', final_x);
fprintf('y = %.1f m\n', final_y);

fprintf('\nOptimal angle:\n');
fprintf('%.1f degrees\n', best_angle);

fprintf('\nMaximum weight:\n');
fprintf('%.6f\n', max_weight);
%%
% plot trajectory
plot_line_trajectory(distance, best_angle, final_x, final_y)



%% plots for the weights
plot_weights_grid(search_history, 4)
