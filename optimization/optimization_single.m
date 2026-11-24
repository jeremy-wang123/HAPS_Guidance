function optimization_single(maxtime, B, Kmax)
%% Optimization Formulation
% maxtime: maximum solver time [seconds]
% B: maximum flight distance [m]
% Kmax: maximum number of flight segments

validateattributes(maxtime, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(B, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(Kmax, {'numeric'}, {'scalar','integer','positive','finite'});

set(0, 'DefaultFigureVisible', 'off');
scriptdir = fileparts(mfilename('fullpath'));
load(fullfile(scriptdir, 'ROI2.mat'));

%% Create flight nodes (V) and evaluation nodes (Q)

g = 5e3;      % Flight node spacing [m]
h = 1e3;      % Evaluation node spacing [m]
buffer = 10e3; %  buffer zone

roi_flight = polybuffer(roi_2, buffer);

%% Radar parameters

% Smaller swath distances for the test problem
R_near = 5e3;      % Near-range distance [m]
R_far  = 10e3;      % Far-range distance [m]

%% Flight nodes (buffered ROI)

% Get buffered ROI boundary
[xb_flight, yb_flight] = boundary(roi_flight);

% Keep the same coordinate origin
[xb, yb] = boundary(roi_2);
x0 = mean(xb);
y0 = mean(yb);

% Transform buffered boundary into strike-normal coordinates
dx_flight = xb_flight - x0;
dy_flight = yb_flight - y0;

s_boundary_flight = dx_flight*strike(1) + dy_flight*strike(2);
n_boundary_flight = dx_flight*normal(1) + dy_flight*normal(2);

% Generate flight nodes
[s_flight, n_flight] = meshgrid( ...
    min(s_boundary_flight):g:max(s_boundary_flight), ...
    min(n_boundary_flight):g:max(n_boundary_flight));

x_flight = x0 + s_flight*strike(1) + n_flight*normal(1);
y_flight = y0 + s_flight*strike(2) + n_flight*normal(2);

% Keep nodes inside buffered ROI
inside_flight = isinterior(roi_flight, x_flight(:), y_flight(:));

V = [x_flight(inside_flight), y_flight(inside_flight)];
N = size(V,1);


%% Evaluation nodes (original ROI)

dx = xb - x0;
dy = yb - y0;

s_boundary = dx*strike(1) + dy*strike(2);
n_boundary = dx*normal(1) + dy*normal(2);

[s_eval, n_eval] = meshgrid( ...
    min(s_boundary):h:max(s_boundary), ...
    min(n_boundary):h:max(n_boundary));

x_eval = x0 + s_eval*strike(1) + n_eval*normal(1);
y_eval = y0 + s_eval*strike(2) + n_eval*normal(2);

inside_eval = isinterior(roi_2, x_eval(:), y_eval(:));

Q = [x_eval(inside_eval), y_eval(inside_eval)];
J = size(Q,1);

%% Plotting
% Display results
fprintf('Flight nodes: %d\n', N);
fprintf('Evaluation nodes: %d\n', J);
fprintf('Candidate directed segments: %d\n', N*(N-1));

% Plot
figure;
hold on;

plot(roi_2, ...
    'FaceColor', [0.9 0.9 0.9], ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'b', ...
    'LineWidth', 1.5);

plot(roi_flight, ...
    'FaceColor', [0.9 0.9 0.9], ...
    'FaceAlpha', 0.2, ...
    'EdgeColor', 'r', ...
    'LineWidth', 1.5);


scatter(Q(:,1), Q(:,2), 3, 'b', 'filled');
scatter(V(:,1), V(:,2), 35, 'r', 'filled');

axis equal;
grid on;
xlabel('x (m)');
ylabel('y (m)');
legend('ROI', 'Flight ROI','Evaluation nodes','Flight nodes');

%{
% Folder path
save_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jeremy.c.wang.27@dartmouth.edu/.shortcut-targets-by-id/13AU22WqA4BTuY5cXBMOREjCRg0XNpKkM/HAPS_Guidance/Figures/optimization_figures';

% Create folder if it doesn't exist
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% Save current figure as PNG (300 dpi)
exportgraphics(gcf, fullfile(save_dir, 'optimization_roi.png'), ...
    'Resolution', 300);
%}
%% Candidate flights
% Generate all directed node pairs
[a,b] = ndgrid(1:N,1:N);

% Remove self-connections
valid = a ~= b;

E = [a(valid), b(valid)];

% Compute segment geometry
p1 = V(E(:,1),:);
p2 = V(E(:,2),:);

d = p2-p1;
flight_length = hypot(d(:,1),d(:,2));

% Maximum allowed individual segment length
Lmax = 30e3;  % [m], adjust based on flight requirements

% Keep only segments within length limit
keep = flight_length <= Lmax;

E = E(keep,:);
p1 = p1(keep,:);
p2 = p2(keep,:);
d = d(keep,:);
flight_length = flight_length(keep);

% Updated number of segments
M = size(E,1);

% Unit flight direction
u = d ./ flight_length;

% Right-facing unit normal
v = [u(:,2), -u(:,1)];

fprintf('Candidate segments after filtering: %d\n',M);

%% Calculate coverage
% A matrix has rows representing each evaluation node, and then columns
% representing each possible candidate flight. It has a binary value for
% whether a candidate flight covers a given node. This is all stored in
% sparse matrix form for efficient memory allocation

batchSize = 100;
numBatches = ceil(M/batchSize);

% Store sparse coverage matrices for each batch
A_blocks = cell(numBatches,1);

tic;

for k = 1:batchSize:M

    % Flight indices for current batch
    idx = k:min(k+batchSize-1,M);

    % Relative positions of evaluation nodes
    dx = Q(:,1) - p1(idx,1)';
    dy = Q(:,2) - p1(idx,2)';

    % Along-track distances
    s = dx .* u(idx,1)' + dy .* u(idx,2)';

    % Rightward cross-track distances
    n = dx .* v(idx,1)' + dy .* v(idx,2)';

    % Check whether nodes fall inside radar swath
    covered = (s >= 0) & ...
              (s <= flight_length(idx)') & ...
              (n >= R_near) & ...
              (n <= R_far);

    % Store coverage for current batch
    batchID = floor((k-1)/batchSize) + 1;
    A_blocks{batchID} = sparse(covered);

    % Progress indicator
    if mod(batchID,100) == 0
        elapsed = toc;
        fraction = idx(end)/M;
        remaining = elapsed/fraction - elapsed;

        fprintf('%.1f%% complete | %.1f min elapsed | %.1f min remaining\n', ...
            100*fraction, elapsed/60, remaining/60);
    end

end

% Combine sparse matrices
A = [A_blocks{:}];

clear A_blocks;

fprintf('Total runtime: %.2f minutes\n',toc/60);
fprintf('Coverage matrix size: %d x %d\n',size(A,1),size(A,2));
fprintf('Nonzero coverage entries: %d\n',nnz(A));
%% Visualize candidate flight coverage

% Choose a candidate segment with the most coverage
coverage_count = full(sum(A,1));
[~,e_test] = max(coverage_count);

% Covered evaluation nodes
covered = full(A(:,e_test)) > 0;

% Flight endpoints
p_start = p1(e_test,:);
p_end = p2(e_test,:);

figure;
hold on;
axis equal;
grid on;

% ROI boundary
plot(roi_2, ...
    'FaceColor',[0.8 0.9 1], ...
    'FaceAlpha',0.1, ...
    'EdgeColor','b');

% Buffered ROI boundary
plot(roi_flight, ...
    'FaceColor','none', ...
    'EdgeColor',[0.6 0.6 0.6], ...
    'LineStyle','--');

% Uncovered evaluation nodes
scatter(Q(~covered,1),Q(~covered,2), ...
    15,[0.7 0.7 0.7],'filled');

% Covered evaluation nodes
scatter(Q(covered,1),Q(covered,2), ...
    35,'g','filled');

% Flight nodes
scatter(V(:,1),V(:,2), ...
    45,'k','filled');

% Selected flight segment
plot([p_start(1),p_end(1)], ...
     [p_start(2),p_end(2)], ...
     'r-','LineWidth',3);

% Flight direction arrow
midpoint = (p_start+p_end)/2;
direction = u(e_test,:);

quiver(midpoint(1),midpoint(2), ...
       direction(1)*3e3,direction(2)*3e3, ...
       0,'r','LineWidth',2, ...
       'MaxHeadSize',2);

% Plot radar swath boundaries
near_start = p_start + R_near*v(e_test,:);
near_end   = p_end   + R_near*v(e_test,:);

far_start  = p_start + R_far*v(e_test,:);
far_end    = p_end   + R_far*v(e_test,:);

% Swath polygon
swathX = [near_start(1),near_end(1), ...
          far_end(1),far_start(1)];

swathY = [near_start(2),near_end(2), ...
          far_end(2),far_start(2)];

patch(swathX,swathY,[0.2 0.8 0.3], ...
    'FaceAlpha',0.15, ...
    'EdgeColor',[0.1 0.6 0.2], ...
    'LineStyle','--');

xlabel('X [m]');
ylabel('Y [m]');

% Include flight endpoints and all four swath corners
x_region = [p_start(1), p_end(1), swathX];
y_region = [p_start(2), p_end(2), swathY];

% Add 5 km padding
padding = 5e3;

xlim([min(x_region)-padding, max(x_region)+padding]);
ylim([min(y_region)-padding, max(y_region)+padding]);

hold off;

%{
% Folder path
save_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jeremy.c.wang.27@dartmouth.edu/.shortcut-targets-by-id/13AU22WqA4BTuY5cXBMOREjCRg0XNpKkM/HAPS_Guidance/Figures/optimization_figures';

% Create folder if it doesn't exist
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% Save current figure as PNG (300 dpi)
exportgraphics(gcf, fullfile(save_dir, 'sample_flight.png'), ...
    'Resolution', 300);
%} 

%% Include all boundary and start and end nodes

%{
% Boundary of buffered ROI
[xb, yb] = boundary(roi_flight);

% Distance from each flight node to boundary
distance_to_boundary = zeros(N,1);

for i = 1:N
    distance_to_boundary(i) = min( ...
        hypot(V(i,1)-xb, V(i,2)-yb));
end

% Boundary tolerance
boundary_tol = g;    % One grid spacing

% Select nodes near boundary, regardless of buffer membership
boundary_nodes = find(distance_to_boundary <= boundary_tol);

% Eligible starting and ending nodes
S = boundary_nodes;
T = boundary_nodes;

%}

%% Include all boundary as start and end nodes 

%{
% Select nodes near the bottom and top of the buffered ROI
start_band = 10e3;  % Bottom 10 km
end_band   = 10e3;  % Top 10 km

% Y coordinates of boundary nodes
y_boundary_nodes = V(boundary_nodes,2);

% Find top and bottom extents
y_min = min(y_boundary_nodes);
y_max = max(y_boundary_nodes);

% Eligible starting nodes: bottom of ROI
S = boundary_nodes(y_boundary_nodes <= y_min + start_band);

% Eligible ending nodes: top of ROI
T = boundary_nodes(y_boundary_nodes >= y_max - end_band);

fprintf('Eligible starting nodes: %d\n',numel(S));
fprintf('Eligible ending nodes: %d\n',numel(T));


% Select nodes near the bottom and top of the buffered ROI
start_band = 10e3;  % Bottom 10 km
end_band   = 10e3;  % Top 10 km

% Y coordinates of boundary nodes
y_boundary_nodes = V(boundary_nodes,2);

% Find top and bottom extents
y_min = min(y_boundary_nodes);
y_max = max(y_boundary_nodes);

% Eligible starting nodes: bottom of ROI
S = boundary_nodes(y_boundary_nodes <= y_min + start_band);

% Eligible ending nodes: top of ROI
T = boundary_nodes(y_boundary_nodes >= y_max - end_band);

fprintf('Eligible starting nodes: %d\n',numel(S));
fprintf('Eligible ending nodes: %d\n',numel(T));

%}
%% Plot start and end nodes

%{
figure;
hold on;
axis equal;
grid on;

plot(roi_flight,'FaceColor','none','EdgeColor','k');
plot(roi_2,'FaceColor','none','EdgeColor','b');

% All flight nodes
scatter(V(:,1),V(:,2),25,[0.7 0.7 0.7],'filled');

% Starting nodes (bottom)
scatter(V(S,1),V(S,2),100,'g','filled');

% Ending nodes (top)
scatter(V(T,1),V(T,2),100,'r','filled');

xlabel('X [m]');
ylabel('Y [m]');
title('Eligible Start and End Nodes');

legend('Buffered ROI','Original ROI', ...
       'Flight Nodes','Starting Nodes','Ending Nodes');

hold off;

%}
%% 

% Boundary of buffered ROI
[xb, yb] = boundary(roi_flight);

% Distance from each flight node to boundary
distance_to_boundary = zeros(N,1);

for i = 1:N
    distance_to_boundary(i) = min( ...
        hypot(V(i,1)-xb, V(i,2)-yb));
end

% Boundary tolerance
boundary_tol = g;    % One grid spacing

% Select nodes near boundary, regardless of buffer membership
boundary_nodes = find(distance_to_boundary <= boundary_tol);

%{
% Select nodes near the bottom and top of the buffered ROI
start_band = 10e3;  % Bottom 10 km
end_band   = 10e3;  % Top 10 km

% Y coordinates of boundary nodes
y_boundary_nodes = V(boundary_nodes,2);

% Find top and bottom extents
y_min = min(y_boundary_nodes);
y_max = max(y_boundary_nodes);

% Eligible starting nodes: bottom of ROI
S = boundary_nodes(y_boundary_nodes <= y_min + start_band);

% Eligible ending nodes: top of ROI
T = S; 
%}

% Select nodes near the bottom of the buffered ROI
start_band = 10e3;

y_boundary_nodes = V(boundary_nodes,2);
y_min = min(y_boundary_nodes);

% Candidate starting nodes
S_candidates = boundary_nodes(y_boundary_nodes <= y_min + start_band);

% Select the bottommost node
[~, idx] = min(V(S_candidates,2));
S = S_candidates(idx);

% Return to the same node
T = S;

fprintf('Eligible starting nodes: %d\n',numel(S));
fprintf('Eligible ending nodes: %d\n',numel(T));

% Plot start and end nodes

figure;
hold on;
axis equal;
grid on;

plot(roi_flight,'FaceColor','none','EdgeColor','k');
plot(roi_2,'FaceColor','none','EdgeColor','b');

% All flight nodes
scatter(V(:,1),V(:,2),25,[0.7 0.7 0.7],'filled');

% Starting nodes (bottom)
scatter(V(S,1),V(S,2),100,'g','filled');

% Ending nodes (top)
scatter(V(T,1),V(T,2),50,'r','filled');

xlabel('X [m]');
ylabel('Y [m]');
title('Eligible Start and End Nodes');

legend('Buffered ROI','Original ROI', ...
       'Flight Nodes','Starting Nodes','Ending Nodes');

hold off;

%% Formulation of optimization constraints
%{
B = 400e3;       % Maximum flight distance [m]
Kmax = 25;       % Maximum number of flight segments

figdir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jeremy.c.wang.27@dartmouth.edu/.shortcut-targets-by-id/13AU22WqA4BTuY5cXBMOREjCRg0XNpKkM/HAPS_Guidance/Figures/optimization_figures/MILP'; 
figname = sprintf('%dB%dK_RTS', B/1e3, Kmax);

[sol, fval, exitflag, output, elapsed_time, ...
 selected_segments, selected_edges, total_distance] = ...
    optimize_flight_MILP(V, E, Q, A, flight_length, ...
                         S, T, B, Kmax, ...
                         roi_2, v, R_near, R_far, true, true, 60*60, ...
                         figdir, figname);
%}


%% Run a single optimization

figdir = fullfile(fileparts(mfilename('fullpath')), 'results_Kmax');

if ~exist(figdir, 'dir')
    mkdir(figdir);
end

    figname = sprintf('%dB%dK_RTS', B/1e3, Kmax);

    % Start recording Command Window output
    logname = fullfile(figdir, [figname, '.txt']);
    diary off;
    diary(logname);
    diary on;

    fprintf('\n\n');
    fprintf('============================================================\n');
    fprintf('============================================================\n');
    fprintf('                 STARTING MILP OPTIMIZATION\n');
    fprintf('                   B = %.0f km | Kmax = %d\n', B/1e3, Kmax);
    fprintf('============================================================\n');
    fprintf('============================================================\n\n');

    [sol, fval, exitflag, output, elapsed_time, ...
     selected_segments, selected_edges, total_distance] = ...
        optimize_flight_MILP(V, E, Q, A, flight_length, ...
                             S, T, B, Kmax, ...
                             roi_2, v, R_near, R_far, true, true, maxtime, ...
                             figdir, figname);

    % Save numerical results
    save(fullfile(figdir, [figname, '.mat']), ...
        'sol', 'fval', 'exitflag', 'output', ...
        'elapsed_time', 'selected_segments', ...
        'selected_edges', 'total_distance', ...
        'B', 'Kmax', 'maxtime', '-v7.3');

    % Stop recording
    diary off;

%% Theoretical parameters

groundspeed = 20; % m/s
max_revisit = 12; % hr
max_distance = groundspeed * 60 * 60 * max_revisit; % m

end % optimization_single
