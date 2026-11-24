%% Test Mixed Integer Linear Program

%% Establish test region of interest

%% Create circular ROI


a = 50e3;              % Semi-major axis [m]
b = 20e3;              % Semi-minor axis [m]
center = [0, 0];

theta = linspace(0, 2*pi, 300);

x_ellipse = center(1) + a*cos(theta);
y_ellipse = center(2) + b*sin(theta);

roi_2 = polyshape(x_ellipse, y_ellipse);


%% Create flight nodes (V) and evaluation nodes (Q)

g = 10e3;              % Flight node spacing [m]
h = 5e3;               % Evaluation node spacing [m]
buffer = 5e3;          % Buffer zone [m]

roi_flight = polybuffer(roi_2, buffer);


%% Flight nodes (buffered ROI)

% Get buffered ROI boundary
[xb_flight, yb_flight] = boundary(roi_flight);

% Keep the same coordinate origin
[xb, yb] = boundary(roi_2);

x0 = center(1);
y0 = center(2);

R_flight = max(hypot(xb_flight-x0, yb_flight-y0));

% Generate grid centered at (x0,y0)
n = ceil(R_flight/g);

[x_flight, y_flight] = meshgrid( ...
    x0 + (-n:n)*g, ...
    y0 + (-n:n)*g);

% Keep nodes inside buffered ROI
inside_flight = isinterior(roi_flight, ...
    x_flight(:), y_flight(:));

V = [x_flight(inside_flight), ...
     y_flight(inside_flight)];

N = size(V,1);


%% Evaluation nodes (original ROI)

% Generate grid covering the original ellipse
[x_eval, y_eval] = meshgrid( ...
    x0-a:h:x0+a, ...
    y0-b:h:y0+b);

% Keep nodes inside original ROI
inside_eval = isinterior(roi_2, ...
    x_eval(:), y_eval(:));

Q = [x_eval(inside_eval), ...
     y_eval(inside_eval)];

J = size(Q,1);

fprintf('Number of evaluation nodes: %d\n', J);

%% Plot ROI and nodes

figure;
hold on;
axis equal;
grid on;

% Buffered ROI
plot(roi_flight, ...
    'FaceColor',[0.9 0.9 0.9], ...
    'FaceAlpha',0.2, ...
    'EdgeColor',[0.5 0.5 0.5]);

% Original ROI
plot(roi_2, ...
    'FaceColor',[0.6 0.8 1], ...
    'FaceAlpha',0.15, ...
    'EdgeColor','b');

% Evaluation nodes
scatter(Q(:,1),Q(:,2), ...
    40,'k','filled');

% Flight nodes
scatter(V(:,1),V(:,2), ...
    100,'r','filled');

xlabel('X [m]');
ylabel('Y [m]');
title('Circular ROI with Flight and Evaluation Nodes');

legend('Buffered ROI','Original ROI', ...
    'Evaluation Nodes','Flight Nodes', ...
    'Location','bestoutside');

hold off;

fprintf('Number of flight nodes: %d\n',N);
fprintf('Number of evaluation nodes: %d\n',J);


%% Candidate flights

% Every directed pair of distinct flight nodes
[a,b] = ndgrid(1:N,1:N);

valid = a ~= b;

E = [a(valid), b(valid)];
M = size(E,1);

%% Radar parameters

% Smaller swath distances for the test problem
R_near = 2e3;      % Near-range distance [m]
R_far  = 8e3;      % Far-range distance [m]

%% Precompute flight geometry

% Starting and ending coordinates
p1 = V(E(:,1),:);
p2 = V(E(:,2),:);

% Flight direction and length
d = p2 - p1;
flight_length = hypot(d(:,1),d(:,2));

% Unit flight direction
u = d ./ flight_length;

% Rightward unit normal
v = [u(:,2), -u(:,1)];

%% Calculate coverage matrix A

batchSize = 100;
numBatches = ceil(M/batchSize);

A_blocks = cell(numBatches,1);

tic;

for k = 1:batchSize:M

    idx = k:min(k+batchSize-1,M);

    % Relative positions of evaluation nodes
    dx = Q(:,1) - p1(idx,1)';
    dy = Q(:,2) - p1(idx,2)';

    % Along-track distances
    s = dx .* u(idx,1)' + ...
        dy .* u(idx,2)';

    % Rightward cross-track distances
    n = dx .* v(idx,1)' + ...
        dy .* v(idx,2)';

    % Radar coverage condition
    covered = (s >= 0) & ...
              (s <= flight_length(idx)') & ...
              (n >= R_near) & ...
              (n <= R_far);

    % Store sparse coverage matrix
    batchID = floor((k-1)/batchSize) + 1;
    A_blocks{batchID} = sparse(covered);

end

% Combine sparse matrices
A = [A_blocks{:}];

clear A_blocks;

fprintf('Total runtime: %.2f seconds\n',toc);
fprintf('Flight nodes: %d\n',N);
fprintf('Evaluation nodes: %d\n',J);
fprintf('Candidate flights: %d\n',M);
fprintf('Coverage matrix size: %d x %d\n',size(A,1),size(A,2));
fprintf('Nonzero entries: %d\n',nnz(A));

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

title(sprintf( ...
    'Candidate Flight %d: %d Evaluation Nodes Covered', ...
    e_test,nnz(covered)));

hold off;

%% Establish boundary nodes

% Define eligible start and end nodes

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

figure;
hold on;
axis equal;
grid on;

plot(roi_flight, ...
    'FaceColor',[0.9 0.9 0.9], ...
    'FaceAlpha',0.1);

plot(roi_2, ...
    'FaceColor',[0.6 0.8 1], ...
    'FaceAlpha',0.1);

% All flight nodes
scatter(V(:,1),V(:,2),40,[0.6 0.6 0.6],'filled');

% Eligible starting/ending nodes
scatter(V(boundary_nodes,1),V(boundary_nodes,2), ...
    100,'r','filled');

xlabel('X [m]');
ylabel('Y [m]');
title('Eligible Start and End Nodes on ROI Boundary');

legend('Buffered ROI','Original ROI', ...
       'Flight Nodes','Boundary Nodes');

hold off;

%% MILP Optimization

B = 80e3;       % Maximum flight distance [m]
Kmax = 10;       % Maximum number of flight segments

[sol, fval, exitflag, output, elapsed_time, ...
 selected_segments, selected_edges, total_distance] = ...
    optimize_flight_MILP(V, E, Q, A, flight_length, ...
                         S, T, B, Kmax, ...
                         roi_2, v, R_near, R_far, true);
