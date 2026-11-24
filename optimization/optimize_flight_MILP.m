
function [sol, fval, exitflag, output, elapsed_time, ...
          selected_segments, selected_edges, total_distance] = ...
          optimize_flight_MILP(V, E, Q, A, flight_length, ...
                               S, T, B, Kmax, ...
                               roi_2, v, R_near, R_far, ...
                               plot_results, return_to_start, maxtime,...
                               figdir, figname)
% OPTIMIZE_FLIGHT_MILP Optimize a flight trajectory for radar coverage.
%
% Solves a Mixed-Integer Linear Program (MILP) that maximizes the
% number of covered evaluation nodes subject to flight continuity,
% starting/ending node restrictions, maximum flight distance, and
% maximum number of flight segments.
%
% INPUTS:
%   V               - [N x 2] Flight node coordinates [m].
%   E               - [M x 2] Directed flight segments. Each row
%                     contains [start_node, end_node] indices into V.
%   Q               - [J x 2] Evaluation node coordinates [m].
%   A               - [J x M] Coverage matrix. A(j,e) = 1 if
%                     segment e covers evaluation node j.
%   flight_length   - [M x 1] Length of each flight segment [m].
%   S               - Vector of eligible starting node indices.
%   T               - Vector of eligible ending node indices.
%   B               - Maximum total flight distance [m].
%   Kmax            - Maximum number of flight segments.
%   roi_2           - polyshape defining the region of interest.
%   v               - [M x 2] Right-facing unit normal vector
%                     for each directed flight segment.
%   R_near          - Near radar swath offset [m].
%   R_far           - Far radar swath offset [m].
%   plot_results    - Logical flag to generate trajectory plot.
%   return_to_start - enforce starting and ending node are the same
%   maxtime         - max time allotted in seconds 
%   figdir          - figure directory
%   figname         - figure name
% 
% OUTPUTS:
%   sol               - Structure containing optimized decision
%                       variables x, y, r, s, t, u_seg, and z.
%   fval              - Optimized objective value.
%   exitflag          - MATLAB solver exit flag.
%   output            - Solver diagnostic information.
%   elapsed_time      - Total optimization runtime [s].
%   selected_segments - [Kused x 1] Selected segment indices.
%   selected_edges    - [Kused x 2] Start/end node indices
%                       of the optimized trajectory.
%   total_distance    - Total optimized flight distance [m].
%
% NOTES:
%   - Requires MATLAB Optimization Toolbox.
%   - The objective uses equal weights for all evaluation nodes.
%   - The MILP is solved using intlinprog.
%   - If no feasible solution is returned, trajectory outputs
%     remain empty and total_distance is NaN.
%
% EXAMPLE:
%   [sol, fval, exitflag, output, elapsed_time, ...
%    selected_segments, selected_edges, total_distance] = ...
%       optimize_flight_MILP(V, E, Q, A, flight_length, ...
%                            S, T, 150e3, 10, ...
%                            roi_2, v, R_near, R_far, true);

%% Initialization

% Optional plotting
if nargin < 14 || isempty(plot_results)
    plot_results = true;
end

% Optional return-to-start constraint
if nargin < 15 || isempty(return_to_start)
    return_to_start = false;
end

if return_to_start && isempty(intersect(S,T))
    error(['Return-to-start requires at least one node ', ...
           'belonging to both S and T.']);
end

N = size(V,1);
M = size(E,1);
J = size(Q,1);


flight_length = flight_length(:);
S = S(:);
T = T(:);

% Initialize outputs
sol = [];
fval = NaN;
exitflag = NaN;
output = [];
elapsed_time = NaN;
selected_segments = [];
selected_edges = [];
total_distance = NaN;


%% Decision Variables

% Flight segment selection
x = optimvar('x',Kmax,M,'Type','integer', ...
    'LowerBound',0,'UpperBound',1);

% Active flight steps
y = optimvar('y',Kmax,1,'Type','integer', ...
    'LowerBound',0,'UpperBound',1);

% Termination variable
r = optimvar('r',Kmax,N,'Type','integer', ...
    'LowerBound',0,'UpperBound',1);

% Starting and ending nodes
s = optimvar('s',N,1,'Type','integer', ...
    'LowerBound',0,'UpperBound',1);

t = optimvar('t',N,1,'Type','integer', ...
    'LowerBound',0,'UpperBound',1);

% Segment utilization
u_seg = optimvar('u_seg',M,1, ...
    'LowerBound',0,'UpperBound',1);

% Evaluation node coverage
z = optimvar('z',J,1, ...
    'LowerBound',0,'UpperBound',1);

%% Section 2: Objective Function

w = ones(J,1);
prob = optimproblem('ObjectiveSense','maximize');
prob.Objective = w' * z;

%% Section 3.1: Max Flight Distance 

% equation 3.1
prob.Constraints.max_distance = ...
    sum(x * flight_length) <= B;

%% Section 3.2: Flight Steps Permitted

% Equation (3.3): One segment per active step
prob.Constraints.one_segment_per_step = ...
    sum(x,2) == y;

% Equation (3.4): Consecutive active steps
if Kmax > 1
    prob.Constraints.consecutive_steps = ...
        y(2:Kmax) <= y(1:Kmax-1);
end

% Equation (3.5): First step must be active
prob.Constraints.first_step = y(1) == 1;

%% Section 3.3: Flight Continuity

% Equation (3.7): Incidence matrices
C_in = sparse(E(:,2), (1:M)', 1, N, M);
C_out = sparse(E(:,1), (1:M)', 1, N, M);

% Equation (3.7): Flight continuity
for k = 1:Kmax-1
    prob.Constraints.(sprintf('continuity_%d',k)) = ...
        C_in*x(k,:)' == ...
        C_out*x(k+1,:)' + r(k,:)';
end

% Equation (3.7): Termination only at T
notT = setdiff(1:N,T);

prob.Constraints.invalid_termination = ...
    r(:,notT) == 0;

% Equation (3.8): Flight termination
prob.Constraints.termination = ...
    sum(r(1:Kmax-1,T),2) == ...
    y(1:Kmax-1) - y(2:Kmax);

%% 3.4: Starting and Ending Points

% Equation (3.9): Exactly one start and end node
prob.Constraints.one_start = sum(s(S)) == 1;
prob.Constraints.one_end = sum(t(T)) == 1;

% Restrict starting/ending nodes to S/T
notS = setdiff(1:N,S);

prob.Constraints.invalid_start = s(notS) == 0;
prob.Constraints.invalid_end = t(notT) == 0;

% Equation (3.10): Ending node selection
prob.Constraints.end_selection = ...
    t(T) == sum(r(:,T),1)';

% Equation (3.11): Starting node selection
prob.Constraints.start_selection = ...
    C_out*x(1,:)' == s;

% Equation (3.12): Final-step termination
prob.Constraints.final_termination = ...
    r(Kmax,:)' == C_in*x(Kmax,:)';

%% 3.4.1: Return to starting node

 % Equation 3.13
if return_to_start
    prob.Constraints.return_to_start = s == t;
end

%% Section 3.5: Flight Segment Selector

% Equation (3.15): Segment utilization
prob.Constraints.segment_usage = ...
    u_seg <= sum(x,1)';

%% Section 3.6: Coverage Constraint

% Equation (3.16): Evaluation node coverage
prob.Constraints.coverage = z <= A * u_seg;

%% Solve MILP

options = optimoptions('intlinprog', ...
    'Display','iter', ...
    'RelativeGapTolerance',0.01, ...
    'MaxTime',maxtime);

tic;

[sol, fval, exitflag, output] = solve(prob, ...
    'Solver','intlinprog', ...
    'Options',options);

elapsed_time = toc;


%% Solver Diagnostics

fprintf('\n--- MILP Optimization Results ---\n');
fprintf('Exit flag: %d\n',exitflag);
fprintf('Solver runtime: %.2f seconds\n',elapsed_time);

if isfield(output,'message')
    disp(output.message);
end

% Stop if the solver returned no usable solution
if isempty(sol) || ~isfield(sol,'x') || ...
        isempty(sol.x) || ~isfield(sol,'y') || ...
        isempty(sol.y)

    fprintf('No feasible trajectory available.\n');
    return;
end

fprintf('Objective value: %.4f\n',fval);
fprintf('Covered evaluation nodes: %d / %d\n', ...
    nnz(sol.z > 0.5),J);


%% Extract Optimized Trajectory

% Equation (1.18): Number of active flight steps
Kused = round(sum(sol.y));

selected_segments = zeros(Kused,1);

for k = 1:Kused

    e = find(sol.x(k,:) > 0.5,1);

    if isempty(e)
        warning('No segment selected at active step %d.',k);
        selected_segments = [];
        selected_edges = [];
        return;
    end

    selected_segments(k) = e;
end

% Selected starting/ending nodes
selected_edges = E(selected_segments,:);

% Total distance
total_distance = sum(flight_length(selected_segments));

fprintf('Flight segments used: %d\n',Kused);
fprintf('Total flight distance: %.2f km\n', ...
    total_distance/1000);

for k = 1:Kused
    e = selected_segments(k);

    start_node = E(e,1);
    end_node = E(e,2);
    segment_length = flight_length(e)/1000;

    fprintf('%d\t%d\t\t%d\t\t%.2f\n', ...
        k, start_node, end_node, segment_length);
end

disp('Selected edges:');
disp(selected_edges);


%% Plot Optimized Trajectory

if plot_results

    figure;
    hold on;
    axis equal;
    grid on;

    % Plot ROI
    plot(roi_2, ...
        'FaceColor','none', ...
        'EdgeColor','k');

    % Covered evaluation nodes
    covered = sol.z > 0.5;

    % Uncovered evaluation nodes
    % scatter(Q(~covered,1),Q(~covered,2), ...
    %     15,[0.7 0.7 0.7],'filled');

    % Plot all flight nodes
    scatter(V(:,1),V(:,2),10,'k','filled');

    % Plot radar swaths
    for k = 1:Kused

        e = selected_segments(k);

        i = E(e,1);
        j = E(e,2);

        p_start = V(i,:);
        p_end = V(j,:);

        % Radar swath boundaries
        near_start = p_start + R_near*v(e,:);
        near_end   = p_end   + R_near*v(e,:);

        far_start = p_start + R_far*v(e,:);
        far_end   = p_end   + R_far*v(e,:);

        % Swath polygon
        swathX = [near_start(1),near_end(1), ...
                  far_end(1),far_start(1)];

        swathY = [near_start(2),near_end(2), ...
                  far_end(2),far_start(2)];

        patch(swathX,swathY,[0.2 0.8 0.3], ...
            'FaceAlpha',0.15, ...
            'EdgeColor',[0.1 0.6 0.2], ...
            'LineStyle','--');
    end

    % Plot selected flight segments
    for k = 1:Kused

        i = selected_edges(k,1);
        j = selected_edges(k,2);

        plot([V(i,1),V(j,1)], ...
             [V(i,2),V(j,2)], ...
             'r-','LineWidth',2);
    end

    % Plot covered evaluation nodes
    scatter(Q(covered,1),Q(covered,2), ...
        5,'m','filled');

    % Mark starting and ending nodes
    if Kused > 0

        start_node = selected_edges(1,1);
        end_node = selected_edges(end,2);

        scatter(V(start_node,1),V(start_node,2), ...
            100,'g','filled');

        scatter(V(end_node,1),V(end_node,2), ...
            100,'b','filled');
    end

    xlabel('X (m)');
    ylabel('Y (m)');
    title(sprintf( ...
        'Optimal Flight Trajectory: %d/%d Nodes Covered', ...
        nnz(covered),J));

    hold off;

    % Create folder if it doesn't exist
    if ~exist(figdir, 'dir')
        mkdir(figdir);
    end
    
    % Save current figure as PNG (300 dpi)
    exportgraphics(gcf, fullfile(figdir, [figname, '.png']), ...
        'Resolution', 300);

end

end
