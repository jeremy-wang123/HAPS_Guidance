%% Elastic Beam Bending varying Young's modulus
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
[h_shelf_mean,h_shelf_eff,h_front_mean,h_front_eff] = extract_shelf_thickness('Thwaites Glacier');
mean_amp = calculate_mean_amp('Thwaites', 2);

% parameters 
rho = 1025; % kg/m^3
g = 9.81; % m/s^2
h = h_front_eff; % use the mean effective thickness of calving front
w0 = mean_amp; % m
sigma_z = 0.0054;  % m (vertical displacement uncertainty)

% parameters for elastic beam bending
I = h^3/12;
mu = 0.325; % poisson's ratio

% Profile specifications
Nel = 160; % number of elements
L = 20e3; % length of ice shelf 
x = linspace(0,L,Nel);
x_km = x/1000;

%% varying E values
E_array = [1e9, 1.4e9, 1.8e9, 2.7e9, 4.7e9, 8.3e9, 9.33e9]; % GPa
sources_array = { ...
    'Vaughan', ...
    'Marsh', ...
    'Saya and Worster', ...
    'Holdsworth', ...
    'Elgart', ...
    'Hooke', ...
    'Gammon'};

%{
Vaughan 1995: 1 GPa
Marsh et al., 2014: 1.4 GPa
Sayag and Worster 2013: 1.8 GPa
Holdsworth 1969: 2.7 GPa
Elgart et al., 2026: 4.7
Hooke (2009): 8.3 GPa
Gammon et al., 1983: 9.33 GPa
%}


figure;
hold on;
for i=1:length(E_array)
    E = E_array(i); 
    D = (E*I)/(1-mu^2); % flexural rigidity

    % Solving w(x) for elastic profile
    lambda = (rho*g*(1-mu^2)/(4*E*I))^(1/4);
    lf = 1/lambda;
    min_L = 5*pi/(4*lambda); % minimum ice shelf length for elastic estimation to make sense (according to Holdsworth)
    w_elastic = w0 * (1 - (exp(-lambda*x) .* (cos(lambda*x) + sin(lambda*x))));

    % elastic
    plot(x_km, w_elastic, ...
        'LineWidth', 2,...
        'DisplayName', sprintf('%s (E = %.2f GPa)', ...
        sources_array{i}, E/1e9)); 

    % Upper and lower uncertainty bounds around steady state
    upper_bound = w0 + sigma_z;
    lower_bound = w0 - sigma_z;

    % Settling-distance markers
    elastic_settling_distance = calculate_settling_distance(w_elastic, x, w0, sigma_z);
    
    align = 'bottom';
    if mod(i,2) == 0
        align = 'middle';
    end

    xline(elastic_settling_distance, ':', ...
        sprintf('E=%.2f GPa: %.1f km', E/1e9, elastic_settling_distance), ...
        'LabelVerticalAlignment', align, ...
        'LabelHorizontalAlignment', 'left', ...
        'HandleVisibility', 'off');
end 

% Far-field steady-state and bounds
yline(w0, '--', '$w_0$', ...
    'Interpreter', 'latex', ...
    'LabelHorizontalAlignment', 'left',...
    'HandleVisibility', 'off');
yline(upper_bound, ':', ...
    '$\pm\sigma_z$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right',...
    'HandleVisibility', 'off');
yline(lower_bound, ':',...
    'HandleVisibility', 'off');

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Deflection at Peak Tide Thwaites Glacier');

xlim([0, L/1000]);

% Axes and labels

xlabel('Distance from grounding line (km)');
ylabel('Vertical displacement (m)');
title('Deflection at Peak Tide: Thwaites Glacier');

xlim([0, max([x_km(:)])]);

% Set y-limits using every plotted profile

all_y_values = [ ...
    w_elastic(:); ...
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

legend('Location', 'southeast');
box on;

set(gca, ...
    'FontSize', 12, ...
    'LineWidth', 1);


exportgraphics(gcf, fullfile(figure_dir,'variable_E_elastic_thwaites.jpg'), ...
'Resolution',300);
