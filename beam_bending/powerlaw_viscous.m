%% Beam Bending
clear; clc;
%% Parameter values 
%{
w0 = tidal amplitude at the edge of the ice shelf (determined from the
CATS2008 model
mu = Poisson's ratio (homogenous for ice)
E = Young's Modulus (uniform for ice)
I = moment of inertia of the section per unit width (see note)
rho = density of water
g = acceleration from gravity

For our purposes, we can treat I  = ∫ y^2 dy from -h/2 to h/2 = h^3 / 12
where h is the thickness of the ice beam. This is assumes uniform thickness
of a rectangular ice shelf and symmetry about y=0

This approximation works for only long shelves (where L >>
5π/4(lambda)

For ice shelves less than that length, we need to treat them as
elastic-plastic deformation
%}

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending';

% extract key params
[h_shelf_mean,h_shelf_eff,h_front_mean,h_front_eff, h_calving_front] = extract_shelf_thickness('Thwaites Glacier');
mean_amp = calculate_mean_amp('Thwaites', 'top2');

% parameters 
rho = 1025; % kg/m^3
g = 9.81; % m/s^2
h = h_calving_front; % use the mean effective thickness of calving front
w0 = mean_amp; % m
sigma_z = 0.0054;  % m (vertical displacement uncertainty)

% parameters for elastic beam bending
I = h^3/12;
mu = 0.325; % poisson's ratio
E = 9.33e9; % Pa
D = (E*I)/(1-mu^2); % flexural rigidity
lambda = (rho*g*(1-mu^2)/(4*E*I))^(1/4);
lf = 1/lambda; % calculated flexural length

% parameters for viscous beam bending
T_hours = 12;            % Tidal period (hours), M2 example
T = T_hours*3600;           % Tidal period (s)
omega = 2*pi/T;             % Angular frequency (rad/s)

A_3 = 3.5e-25; % Pa^-3 s^-1
tau_e_ref = 0.1e6;  
eta_newt = 1/(2*A_3*tau_e_ref^2);
D_v = eta_newt*h^3/3;

% parameters for non-Newtonian FEM
R_gas_local = 8.314; % Gas constant
eps_bg = 1e-10; % background extensional strain rate, 1/s
kappa_bg = eps_bg/(h/2); % equivalent background curvature rate
T_ice = 263.0;  

% flow law parameters
A_1 = 1/(2*eta_newt); % newtonian flow parameter

% n = 4: dislocation creep
A_4 = A_Glen_from_GK( ...
    4.0e5, ...       % A0 in MPa^-n s^-1
    4.0, ...         % n
    60e3, ...        % Q, J/mol
    T_ice, ...
    [], ...
    0, ...
    R_gas_local);

% n = 1.8: GBS-accommodated basal slip
A_18 = A_Glen_from_GK( ...
    3.9e-3, ...      % A0
    1.8, ...         % n
    49e3, ...        % Q, J/mol
    T_ice, ...
    1e-3, ...        % Grain size, m
    1.4, ...         % Grain-size exponent
    R_gas_local);


%% Profile specifications
n_lf = 20;

Nel = 160; % number of elements
L = n_lf*lf; % length of ice shelf 
x = linspace(0,L,Nel);
x_km = x/1000;

%% Solving w(x) for elastic profile

min_L = 5*pi/(4*lambda); % minimum ice shelf length for elastic estimation to make sense (according to Holdsworth)
w_elastic = w0 * (1 - (exp(-lambda*x) .* (cos(lambda*x) + sin(lambda*x))));
%% solving viscous profile

% calculate beta 
beta = (rho*g/(omega*D_v))^(1/4);

% Characteristic viscous flexural length
ell_v = 1/beta;

% characteristic roots
lambda_all = beta*exp(1i*(pi/8 + (0:3)*pi/2));
lambda_decay = lambda_all(real(lambda_all) < 0);
lambda1 = lambda_decay(1);
lambda2 = lambda_decay(2);

% calculates constants
C1 =  w0*lambda2/(lambda1-lambda2);
C2 = -w0*lambda1/(lambda1-lambda2);

% complex spatial response
W = w0 ...
    + C1.*exp(lambda1.*x) ...
    + C2.*exp(lambda2.*x);

% extract local amplitude and phase
amplitude = abs(W);                 % Local tidal amplitude (m)

% W = amplitude*exp(-i*phase_lag)
phase_lag = -angle(W);              % Phase lag (rad)
phase_lag = unwrap(phase_lag);      % Remove artificial 2*pi jumps
phase_lag_deg = rad2deg(phase_lag); % Phase lag (degrees)

% The phase is not meaningful exactly where amplitude is zero.
phase_lag_deg(amplitude < 1e-8*w0) = NaN;

% Peak high tide occurs at t = T/4
t_peak = T/4;

% Instantaneous physical displacement profile
w_peak = imag(W .* exp(1i*omega*t_peak));

%% Power Law viscous tidal flexure model
% ================================================================
%  Flow-law cases
% ================================================================
n_values = [1.0, 1.8, 3.0, 4.0];
A_values = [A_1, A_18, A_3, A_4];

flow_labels = { ...
    'n = 1 (Newtonian)', ...
    'n = 1.8 (GBS)', ...
    'n = 3 (Glen)', ...
    'n = 4 (dislocation)'};

for j = 1:numel(n_values)

    n_val = n_values(j);
    A_val = A_values(j);

    eps_dot_ref = A_val*tau_e_ref^n_val;

    eta_eff_ref = ...
        1/(2*A_val*tau_e_ref^(n_val - 1));

    fprintf(['  %s: A = %.2e, ' ...
             'eps_dot(0.1 MPa) = %.2e 1/s, ' ...
             'eta = %.2e Pa s\n'], ...
        flow_labels{j}, A_val, eps_dot_ref, eta_eff_ref);

end

% ================================================================
% Spatial grid
% ================================================================

he = L/Nel;        % Element length, m

x_nodes = linspace(0, L, Nel + 1)';

% Each node has:
%   1. vertical velocity
%   2. slope rate
ndof = 2*(Nel + 1);


% ================================================================
% Two-point Gauss quadrature
% ================================================================

gp = [ ...
    0.5 - sqrt(3)/6, ...
    0.5 + sqrt(3)/6];

gw = [0.5, 0.5];


% Hermite shape functions

N_gp = zeros(2,4);
d2N_gp = zeros(2,4);

for q = 1:2

    xi = gp(q);

    % Cubic Hermite shape functions
    N_gp(q,:) = [ ...
        1 - 3*xi^2 + 2*xi^3, ...
        xi*(1 - xi)^2, ...
        3*xi^2 - 2*xi^3, ...
        xi^2*(xi - 1)];

    % Second derivatives with respect to xi
    d2N_gp(q,:) = [ ...
        -6 + 12*xi, ...
        -4 + 6*xi, ...
         6 - 12*xi, ...
        -2 + 6*xi];

end


% Convert shape functions to physical coordinates

% The slope degrees of freedom require factors of element length.
scale = [1, he, 1, he];

N_phys = N_gp .* scale;

% d^2/dx^2 = (1/he^2) d^2/dxi^2
B_phys = (d2N_gp/he^2) .* scale;


% ================================================================
% Reference element matrices
% ================================================================

% K_ref(:,:,q) contains:
%
%     weight * he * B' * B
%
% for Gauss point q.

K_ref = zeros(4,4,2);

% F_ref(:,q) contains:
%
%     weight * he * N'
%
F_ref = zeros(4,2);

% Hydrostatic restoring matrix:
%
%     integral N' N dx

M_ref = zeros(4,4);

for q = 1:2

    Bq = B_phys(q,:);
    Nq = N_phys(q,:);

    K_ref(:,:,q) = ...
        gw(q)*he*(Bq'*Bq);

    F_ref(:,q) = ...
        gw(q)*he*Nq';

    M_ref = M_ref + ...
        gw(q)*he*(Nq'*Nq);

end


% ================================================================
% Element-to-global degree-of-freedom mapping
% ================================================================

element_dofs = zeros(Nel,4);

% These arrays are used for efficient sparse-matrix assembly.
row_indices = zeros(Nel,4,4);
col_indices = zeros(Nel,4,4);

for e = 1:Nel

    % Element e connects node e and node e+1.
    %
    % MATLAB degree-of-freedom ordering:
    %
    % node 1: [w_dot_1, slope_dot_1]
    % node 2: [w_dot_2, slope_dot_2]
    %
    dofs = [ ...
        2*e - 1, ...
        2*e, ...
        2*e + 1, ...
        2*e + 2];

    element_dofs(e,:) = dofs;

    row_indices(e,:,:) = repmat(dofs',1,4);
    col_indices(e,:,:) = repmat(dofs,4,1);

end

row_indices_flat = row_indices(:);
col_indices_flat = col_indices(:);

% Store all model information in one structure
model.h = h;
model.w0 = w0;
model.T = T;
model.omega = omega;

model.rho = rho;
model.g = g;

model.Nel = Nel;
model.he = he;
model.ndof = ndof;
model.x_nodes = x_nodes;

model.gp = gp;
model.gw = gw;

model.N_phys = N_phys;
model.B_phys = B_phys;

model.K_ref = K_ref;
model.F_ref = F_ref;
model.M_ref = M_ref;

model.element_dofs = element_dofs;
model.row_indices_flat = row_indices_flat;
model.col_indices_flat = col_indices_flat;

model.kappa_bg = kappa_bg;
model.D_v = D_v;


% ================================================================
% Version 1: Physical flow-law values
% ================================================================

fprintf('\n--- Version 1: Physical A values ---\n');

nl_phys = struct();

for j = 1:numel(n_values)

    n_val = n_values(j);
    A_val = A_values(j);

    D_n = compute_D_n(A_val, n_val, h);

    result = run_power_law( ...
        n_val, ...
        D_n, ...
        flow_labels{j}, ...
        6, ...
        model);

    field_name = matlab.lang.makeValidName( ...
        sprintf('n_%g',n_val));

    nl_phys.(field_name) = result;

end


% ================================================================
% Version 2: Match D_eff(kappa_bg) to D_v
% ================================================================

fprintf('\n--- Version 2: Matched D_eff(kappa_bg) = D_v ---\n');

nl_matched = struct();

for j = 1:numel(n_values)

    n_val = n_values(j);

    % D_eff = D_n*kappa_bg^(1/n - 1)
    %
    % Setting D_eff(kappa_bg) = D_v gives:
    %
    % D_n = D_v*kappa_bg^(1 - 1/n)

    D_n_matched = ...
        D_v*kappa_bg^(1 - 1/n_val);

    matched_label = ...
        sprintf('%s, matched',flow_labels{j});

    result = run_power_law( ...
        n_val, ...
        D_n_matched, ...
        matched_label, ...
        6, ...
        model);

    field_name = matlab.lang.makeValidName( ...
        sprintf('n_%g',n_val));

    nl_matched.(field_name) = result;

end

fprintf('\nPower-law viscous simulations complete.\n');

% Extract results
x_powerlaw_km = x_nodes/1000;

w_n1 = nl_phys.n_1.high_tide_profile;
w_n18 = nl_phys.n_1_8.high_tide_profile;
w_n3  = nl_phys.n_3.high_tide_profile;
w_n4  = nl_phys.n_4.high_tide_profile;

% Make sure all profiles are column vectors
w_n1 = w_n1(:); 
w_n18 = w_n18(:);
w_n3  = w_n3(:);
w_n4  = w_n4(:);


%% Settling-distance calculations

% calculate settling distance analytically
elastic_settling_distance = calculate_settling_distance(w_elastic, x, w0, sigma_z);

% calculate settling distance
viscous_settling_distance = calculate_settling_distance(w_peak, x, w0, sigma_z);

% Power-law FEM profiles
settling_distance_n1 = ...
    calculate_settling_distance( ...
        w_n1, x_nodes, w0, sigma_z);

settling_distance_n18 = ...
    calculate_settling_distance( ...
        w_n18, x_nodes, w0, sigma_z);

settling_distance_n3 = ...
    calculate_settling_distance( ...
        w_n3, x_nodes, w0, sigma_z);

settling_distance_n4 = ...
    calculate_settling_distance( ...
        w_n4, x_nodes, w0, sigma_z);

%% Plotting Profiles

figure;
hold on;
% elastic
plot(x_km, w_elastic, 'k', 'LineWidth', 2, 'DisplayName', 'Elastic'); 

% viscous
plot(x_km, w_peak, ...
    'Color', '#2ca02c', ...
    'LineStyle', '--', ...
    'LineWidth', 2, ...
    'DisplayName', 'Newtonian (analytical)'); 

% Power-law FEM profiles
plot(x_powerlaw_km, w_n18, ...
    'Color', '#d62728',...
    'LineStyle', '-', ...
    'LineWidth', 2, ...
    'DisplayName', '$n=1.8$ power law');

plot(x_powerlaw_km, w_n3, ...
    'Color', '#9467bd',...
    'LineStyle', '-', ...
    'LineWidth', 2, ...
    'DisplayName', '$n=3$ power law');

plot(x_powerlaw_km, w_n4, ...
    'Color', '#1f77b4',...
    'LineStyle', '-', ...
    'LineWidth', 2, ...
    'DisplayName', '$n=4$ power law');

% Far-field steady-state tidal displacement
yline(w0, '--', '$w_0$', ...
    'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left');

% Upper and lower uncertainty bounds around steady state
upper_bound = w0 + sigma_z;
lower_bound = w0 - sigma_z;

yline(upper_bound, ':', ...
    '$\pm\sigma_z$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right');
yline(lower_bound, ':');

% Settling-distance markers
xline(elastic_settling_distance, ':', ...
    sprintf('Elastic = %.1f km', elastic_settling_distance), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left');

xline(viscous_settling_distance, ':', ...
    sprintf('Newtonian = %.1f km', viscous_settling_distance), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left');

xline(settling_distance_n18, ':', ...
    sprintf('n=1.8: %.1f km', settling_distance_n18), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left');

xline(settling_distance_n3, ':', ...
    sprintf('n=3: %.1f km', settling_distance_n3), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left');

xline(settling_distance_n4, ':', ...
    sprintf('n=4: %.1f km', settling_distance_n4), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Deflection at Peak Tide Thwaites Glacier');

xlim([0, L/1000]);

% Axes and labels

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Deflection at Peak Tide: Thwaites Glacier');

xlim([0, max([x_km(:); x_powerlaw_km(:)])]);

% Set y-limits using every plotted profile

all_y_values = [ ...
    w_elastic(:); ...
    w_peak(:); ...
    w_n1(:);...
    w_n18(:); ...
    w_n3(:); ...
    w_n4(:); ...
    upper_bound; ...
    lower_bound];

w_range = max(all_y_values) - min(all_y_values);

if w_range == 0
    vertical_padding = ...
        max(abs(all_y_values))*0.1;
else
    vertical_padding = ...
        0.20*w_range;
end

ylim([ ...
    min(all_y_values) - vertical_padding, ...
    max(all_y_values) + vertical_padding]);


% Plot formatting

legend('Elastic', 'Newtonian', '$n=1.8$ power law', ...
    '$n=3$ power law', '$n=4$ power law',...
    'Location', 'southeast', ...
    'Interpreter', 'latex');

box on;

set(gca, ...
    'FontSize', 12, ...
    'LineWidth', 1);

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/fem';
exportgraphics(gcf, fullfile(figure_dir,'beam_bending_comparison_thwaites.jpg'), ...
'Resolution',300);

%% plotting curvature at peak tide

% analytically solved elastic curvature
elastic_curvature = 2*w0*lambda^2 * exp(-lambda*x) .* (cos(lambda*x) - sin(lambda*x));

% viscous curvature
viscous_slope = gradient(w_peak, x);          % dw/dx
viscous_curvature = gradient(viscous_slope, x);  % d^2w/dx^2

figure;
hold on;
plot(x_km, viscous_curvature, 'LineWidth', 2);
plot(x_km,elastic_curvature, 'LineWidth', 2);

xlabel('Distance from grounding line (km)');
ylabel('Curvature, $d^2w/dx^2 (m^{-1})$', 'Interpreter', 'latex');
title('Curvature at peak tide');

yline(0, 'k--');
xlim([0, L/1000]);

legend('elastic','viscous')
grid on;
box on;
set(gca, 'FontSize', 12, 'LineWidth', 1);

% exportgraphics(gcf, fullfile(figure_dir,'beam_curvature.jpg'), ...
% 'Resolution',300);

%% Deflection at lf

% flexural wave-length (features smaller than lf are support by rigidity,
% features larger than lf are at isostatic equilibrium)

% Find closest grid point
[~,idx] = min(abs(x-lf));

% Time vector (two tidal cycles)
Nt = 500;
t = linspace(0,2*T,Nt);

% Initialize
w_lf = zeros(size(t));

% Displacement at x = lf
for i = 1:Nt

    % Instantaneous displacement profile
    w = imag(W .* exp(1i*omega*t(i)));

    % Record displacement at x = lf
    w_lf(i) = w(idx);

end

% Applied tidal forcing
w_tide = w0*sin(omega*t);

% Plot
figure;
hold on;

plot(t/3600, w_lf, 'LineWidth', 2, ...
    'DisplayName','viscous');

plot(t/3600, w_tide, '--k', 'LineWidth', 2, ...
    'DisplayName','$w_{tide}$');

xlabel('Time (hours)');
ylabel('Vertical displacement (m)');
title(sprintf('Response at $x=\\ell_f = %.2f$ km', lf/1000), ...
    'Interpreter','latex');

legend('Interpreter','latex','Location','best');

grid on;
box on;
set(gca,'FontSize',12,'LineWidth',1);


% figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending';
% exportgraphics(gcf, fullfile(figure_dir,'viscous_phase_lag.jpg'), ...
% 'Resolution',300);

%% Check FEM solver with viscous solution

figure;
hold on;

% viscous
plot(x_km, w_peak, ...
    'Color', '#2ca02c', ...
    'LineStyle', '-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Analytical'); 

% Power-law FEM profiles
plot(x_powerlaw_km, w_n1, ...
    'Color', 'b',...
    'LineStyle', '-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Numerical');

% Far-field steady-state tidal displacement
yline(w0, '--', '$w_0$', ...
    'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

% Upper and lower uncertainty bounds around steady state
upper_bound = w0 + sigma_z;
lower_bound = w0 - sigma_z;

yline(upper_bound, ':', ...
    '$\pm\sigma_z$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right', ...
    'HandleVisibility', 'off');
yline(lower_bound, ':', ...
    'HandleVisibility', 'off');

% Settling-distance markers

xline(viscous_settling_distance, ':', ...
    sprintf('Analytical = %.1f km', viscous_settling_distance), ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xline(settling_distance_n1, ':', ...
    sprintf('Numerical: %.1f km', settling_distance_n1), ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left', ...
    'HandleVisibility', 'off');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Comparison between Numerical and Analytical Viscous Beam');

xlim([0, L/1000]);

% Axes and labels

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Deflection at Peak Tide: Thwaites Glacier');

xlim([0, max([x_km(:); x_powerlaw_km(:)])]);

% Set y-limits using every plotted profile

all_y_values = [ ...
    w_peak(:); ...
    w_n1(:); ...
    upper_bound; ...
    lower_bound];

w_range = max(all_y_values) - min(all_y_values);

if w_range == 0
    vertical_padding = ...
        max(abs(all_y_values))*0.1;
else
    vertical_padding = ...
        0.20*w_range;
end

ylim([ ...
    min(all_y_values) - vertical_padding, ...
    max(all_y_values) + vertical_padding]);


% Plot formatting
legend('Location', 'southeast');

box on;

set(gca, ...
    'FontSize', 12, ...
    'LineWidth', 1);

figure_dir = '/Users/jeremywang/Library/CloudStorage/GoogleDrive-jcwang2@caltech.edu/My Drive/HAPS_Guidance/Figures/beam_bending/viscous_comparisons';
exportgraphics(gcf, fullfile(figure_dir,sprintf('viscous_comparison_%d_lf.jpg', n_lf)), ...
'Resolution',300);
%% Functions for FEM solver

function A_Glen = A_Glen_from_GK( ...
    A0_MPa, n, Q, T, d, p, R_gas)
% A_GLEN_FROM_GK
% Apply temperature, grain-size, and stress-unit conversions.
%
% A0_MPa:
%     Pre-exponential factor with stress measured in MPa.
%
% d:
%     Grain size in meters. Pass [] when grain size is not used.
%
% Important:
% The final GK-to-Glen multiplicative factor should be independently
% verified against the tensor conventions used in the model.

    % Arrhenius temperature dependence
    A_GK_MPa = ...
        A0_MPa*exp(-Q/(R_gas*T));

    % Grain-size dependence
    if ~isempty(d) && p > 0
        A_GK_MPa = A_GK_MPa/d^p;
    end

    % Convert MPa^-n to Pa^-n
    A_GK_Pa = ...
        A_GK_MPa/(1e6)^n;

    % Convention conversion used in the Python code
    A_Glen = ...
        3^((n + 1)/2)/2*A_GK_Pa;

end

function D_n = compute_D_n(A_Glen, n_nl, h_plate)
% COMPUTE_D_N
% Nonlinear viscous bending coefficient for the relation
%
% M = D_n*|kappa_dot|^(1/n - 1)*kappa_dot

    D_n = ...
        (4*n_nl/(1 + 2*n_nl)) ...
        *(h_plate/2)^(1/n_nl + 2) ...
        /A_Glen^(1/n_nl);

end



function wdot_vec = solve_wdot_fem( ...
    w_curr, ...
    w_tide_scalar, ...
    D_n, ...
    n_nl, ...
    dt, ...
    wdot_previous, ...
    picard_iterations, ...
    model)
% SOLVE_WDOT_FEM
% Solve for nodal vertical velocities and nodal slope rates during one
% time step.
%
% Global unknown ordering:
%
% [w_dot_1;
%  slope_dot_1;
%  w_dot_2;
%  slope_dot_2;
%  ... ]

    Nel = model.Nel;
    ndof = model.ndof;

    B_phys = model.B_phys;
    K_ref = model.K_ref;
    F_ref = model.F_ref;
    M_ref = model.M_ref;

    rho = model.rho;
    g = model.g;

    gp = model.gp;

    kappa_bg = model.kappa_bg;

    element_dofs = model.element_dofs;

    row_indices_flat = model.row_indices_flat;
    col_indices_flat = model.col_indices_flat;


    % Initial Picard guess

    if isempty(wdot_previous)
        wdot_vec = zeros(ndof,1);
    else
        wdot_vec = wdot_previous;
    end


    % Linear case requires only one matrix solve

    if n_nl > 1
        number_iterations = picard_iterations;
    else
        number_iterations = 1;
    end


    % Picard nonlinear iteration

    for iteration = 1:number_iterations

        % One effective bending coefficient at each Gauss point
        % of every element.
        D_eff = D_n*ones(Nel,2);


        % Calculate nonlinear effective bending resistance

        if n_nl > 1

            for e = 1:Nel

                dofs = element_dofs(e,:);

                % Element velocity degrees of freedom:
                %
                % [left velocity;
                %  left slope rate;
                %  right velocity;
                %  right slope rate]
                dv = wdot_vec(dofs);

                for q = 1:2

                    % Curvature rate:
                    %
                    % kappa_dot = d^2(w_dot)/dx^2
                    kappa_dot = ...
                        B_phys(q,:)*dv;

                    % Smooth regularization
                    kappa_dot_eff = ...
                        sqrt(kappa_dot^2 + kappa_bg^2);

                    % Nonlinear effective bending coefficient
                    D_eff(e,q) = ...
                        D_n ...
                        *kappa_dot_eff^(1/n_nl - 1);

                end

            end

        end


        % Assemble global matrix element by element

        K_global = sparse(ndof,ndof);
        
        for e = 1:Nel
        
            % Four global degrees of freedom for this element
            dofs = element_dofs(e,:);
        
            % Local 4-by-4 element matrix
            K_e = ...
                  D_eff(e,1)*K_ref(:,:,1) ...
                + D_eff(e,2)*K_ref(:,:,2) ...
                + dt*rho*g*M_ref;
        
            % Add local matrix into the global matrix
            K_global(dofs,dofs) = ...
                K_global(dofs,dofs) + K_e;
        
        end


        % Current displacement at Gauss points

        w_gp = zeros(Nel,2);

        for q = 1:2

            xi = gp(q);

            % This reproduces the Python code exactly:
            % linear interpolation between nodal displacements.
            w_gp(:,q) = ...
                (1 - xi)*w_curr(1:end-1) ...
                + xi*w_curr(2:end);

        end


        % Hydrostatic load

        load_gp = ...
            rho*g*(w_tide_scalar - w_gp);


        % Element force vectors

        F_element = zeros(Nel,4);

        for e = 1:Nel

            F_element(e,:) = ...
                (load_gp(e,1)*F_ref(:,1) ...
                + load_gp(e,2)*F_ref(:,2))';

        end


        % Assemble global force vector

        F_global = accumarray( ...
            element_dofs(:), ...
            F_element(:), ...
            [ndof,1], ...
            @sum, ...
            0);


        % Clamped grounding-line boundary conditions

        constrained_dofs = [1,2];

        % Degree of freedom 1:
        %     w_dot(0) = 0
        %
        % Degree of freedom 2:
        %     d(w_dot)/dx at x=0 = 0

        for bc = constrained_dofs

            K_global(bc,:) = 0;
            K_global(:,bc) = 0;

            K_global(bc,bc) = 1;
            F_global(bc) = 0;

        end


        % Solve the linear system

        wdot_vec = K_global\F_global;

    end

end



function result = run_power_law( ...
    n_nl, ...
    D_n, ...
    label, ...
    number_cycles, ...
    model)
% RUN_POWER_LAW
% Run a power-law viscous simulation for multiple tidal cycles.

    kappa_bg = model.kappa_bg;
    D_v = model.D_v;

    Nel = model.Nel;
    ndof = model.ndof;

    T = model.T;
    omega = model.omega;
    w0 = model.w0;

    he = model.he;


    % Background effective rigidity

    D_eff_background = ...
        D_n*kappa_bg^(1/n_nl - 1);

    fprintf(['  %s: D_n = %.2e, ' ...
             'D_eff(kappa_bg) = %.2e, ' ...
             'ratio to D_v = %.2f\n'], ...
        label, ...
        D_n, ...
        D_eff_background, ...
        D_eff_background/D_v);


    % Time discretization

    dt = 30.0;

    steps_per_cycle = ...
        round(T/dt);

    number_steps = ...
        number_cycles*steps_per_cycle;

    % Save approximately 200 profiles during the last cycle.
    save_every = ...
        max(1,floor(steps_per_cycle/200));

    save_start = ...
        (number_cycles - 1)*steps_per_cycle + 1;


    % Initial conditions

    % Only nodal displacement is explicitly stored here.
    w = zeros(Nel + 1,1);

    wdot_previous = [];

    maximum_saved_profiles = ...
        ceil(steps_per_cycle/save_every) + 1;

    w_history = ...
        zeros(maximum_saved_profiles,Nel + 1);

    time_history = ...
        zeros(maximum_saved_profiles,1);

    save_count = 0;


    % Time stepping

    for step = 1:number_steps

        % MATLAB indexing begins at 1, but physical time begins at zero.
        current_time = ...
            (step - 1)*dt;

        w_tide = ...
            w0*sin(omega*current_time);


        % Solve for vertical velocity

        wdot_vec = solve_wdot_fem( ...
            w, ...
            w_tide, ...
            D_n, ...
            n_nl, ...
            dt, ...
            wdot_previous, ...
            3, ...
            model);


        % Save velocity as next initial Picard guess

        wdot_previous = wdot_vec;


        % Update nodal displacement

        % Odd MATLAB indices 1,3,5,... are vertical velocities.
        vertical_velocity = ...
            wdot_vec(1:2:ndof);

        w = ...
            w + dt*vertical_velocity;

        % Explicitly enforce grounding-line displacement.
        w(1) = 0;


        % Save only the final cycle

        if step >= save_start && ...
                mod(step - save_start,save_every) == 0

            save_count = save_count + 1;

            w_history(save_count,:) = w';
            time_history(save_count) = current_time;

        end

    end


    % Remove unused preallocated rows

    w_history = ...
        w_history(1:save_count,:);

    time_history = ...
        time_history(1:save_count);


    % Find profile nearest positive high tide

    phase = ...
        omega*time_history;

    [~,high_tide_index] = ...
        max(sin(phase));

    high_tide_profile = ...
        w_history(high_tide_index,:);


    % Diagnostics

    maximum_displacement = ...
        max(abs(w_history),[],'all');

    offshore_displacement = ...
        high_tide_profile(end);

    offshore_slope = ...
        (high_tide_profile(end) ...
        - high_tide_profile(end - 1))/he;

    fprintf(['    max|w| = %.3f m, ' ...
             'w(L) = %.4f m, ' ...
             'dw/dx(L) = %.2e\n'], ...
        maximum_displacement, ...
        offshore_displacement, ...
        offshore_slope);

    % Return results

    result.n = n_nl;
    result.D_n = D_n;

    result.w_history = w_history;
    result.time_history = time_history;

    result.final_displacement = w;
    result.final_wdot = wdot_previous;

    result.high_tide_profile = ...
        high_tide_profile';

end