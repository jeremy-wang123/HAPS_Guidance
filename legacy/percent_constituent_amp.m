%% Plot percentage of amplitude top (x) constituents contribute
clear; clc;
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008')
addpath('/Users/jeremywang/Documents/MATLAB/measures_v3.1.2/measures')
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008';

% Thwaites region
% latrange = -78.0:0.02:-73;
% lonrange = -111:0.02:-99;

namps = 2; % how many of the top amps are being taken

latrange = -90.0:0.05:-60;
lonrange = -180:0.05:180;

[lon,lat] = meshgrid(lonrange,latrange);

% Check constituent list
conlist = extract_conlist(Model);

% extracting percentages that top 2 amplitudes constitute
[amp,~,~,~] = tmd_extract_HC(Model,lat,lon,'z');
percentages = NaN(length(latrange), length(lonrange));
for i=1:length(latrange)
    for j=1:length(lonrange)
        amps_per_grid = amp(:,i,j);
        sorted_amps = sort(amps_per_grid);
        topAmps = sorted_amps(end-(namps-1):end); % Get top amplitudes
        totalAmp = sum(amps_per_grid);
        percentages(i,j) = (sum(topAmps, 'omitnan') / totalAmp) * 100; % Calculate percentage contribution
    end
end

%%
% Convert CATS lon/lat grid to polar stereographic x/y used by CATS
[x,y] = ll2ps(lat,lon);

figure
% ice speed plot
ax1 = axes;
% mapzoomps('Thwaites Glacier')
measuresps('speed','log')
measuresps('gl', 'k')
hold on

axis(ax1,'equal')
axis(ax1,'manual')

xl = xlim(ax1);
yl = ylim(ax1);

colormap(ax1, parula)
cb1 = colorbar(ax1,'westoutside');
cb1.Ticks = [0 1 2 3];
cb1.TickLabels = {'1','10','100','1000'};
ylabel(cb1,'Ice speed (m/yr)')

xlabel('x (m)');
ylabel('y (m)');
% percentages plot
ax2 = axes;

pcolor(ax2, x, y, percentages)
shading(ax2,'interp')
set(findobj(ax2,'Type','Surface'), ...
    'FaceAlpha',0.5, ...
    'EdgeColor','none')

colormap(ax2, turbo)
cb2 = colorbar(ax2,'eastoutside');
ylabel(cb2,'Percentage of total amplitude attributed to top 2 constituents')
% caxis([0 100])

% Match axes exactly
set(ax2, ...
    'Position', ax1.Position, ...
    'Color','none', ...
    'XLim', xl, ...
    'YLim', yl, ...
    'DataAspectRatio', ax1.DataAspectRatio, ...
    'PlotBoxAspectRatio', ax1.PlotBoxAspectRatio, ...
    'XTick', [], ...
    'YTick', [], ...
    'Box','off')

linkaxes([ax1 ax2],'xy')

% Put ice-speed axes visually underneath
uistack(ax1,'bottom')

title(ax1,sprintf('Percentage of total amplitude attributed to top %d constituents', namps))

%% Let's try plotting the distribution

percentages_2 = extract_percentages(2);
percentages_3 = extract_percentages(3);
percentages_4 = extract_percentages(4);
percentages_5 = extract_percentages(5);
percentages_6 = extract_percentages(6);
percentages_7 = extract_percentages(7);
percentages_8 = extract_percentages(8);

figure;
histogram(percentages_2, 'Normalization', 'pdf')
hold on; 
histogram(percentages_3, 'Normalization', 'pdf')
histogram(percentages_4, 'Normalization', 'pdf')
histogram(percentages_5, 'Normalization', 'pdf')
histogram(percentages_6, 'Normalization', 'pdf')
histogram(percentages_7, 'Normalization', 'pdf')
histogram(percentages_8, 'Normalization', 'pdf')

xlabel('Percent')
ylabel('Probability Density')
legend('Top 2', 'Top 3', 'Top 4', 'Top 5', 'Top 6', 'Top 7', 'Top 8');
title('Percentage contribution to total amplitude from constituents')

%% function for extracting percentages


function percentages = extract_percentages(namps)
    addpath('/Users/jeremywang/Documents/MATLAB/CATS2008')
    addpath('/Users/jeremywang/Documents/MATLAB/measures_v3.1.2/measures')
    Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008';

    latrange = -90.0:0.05:-60;
    lonrange = -180:0.05:180;
    
    [lon,lat] = meshgrid(lonrange,latrange);
    
    % Check constituent list
    conlist = extract_conlist(Model);
    
    % extracting percentages that top 2 amplitudes constitute
    [amp,~,~,~] = tmd_extract_HC(Model,lat,lon,'z');
    percentages = NaN(length(latrange), length(lonrange));
    for i=1:length(latrange)
        for j=1:length(lonrange)
            amps_per_grid = amp(:,i,j);
            sorted_amps = sort(amps_per_grid);
            topAmps = sorted_amps(end-(namps-1):end); % Get top amplitudes
            totalAmp = sum(amps_per_grid);
            percentages(i,j) = (sum(topAmps, 'omitnan') / totalAmp) * 100; % Calculate percentage contribution
        end
    end
end