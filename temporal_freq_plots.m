%% plotting flight temporal resolution required
% revisit frequency: period of revisit time for each time location based on
% dominant constituent
% campaign time: the amount of time needed to capture the beat frequency
% of the dominant two constituents

%% Plotting dominant constituent
clear;clc;

% loading CATS2008 Model
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008');
addpath('/Users/jeremywang/Documents/MATLAB/measures_v3.1.2/measures');
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008'; 

% mapping region
latrange = -90.0:0.02:-60;
lonrange = -180:0.02:180;
[lon,lat] = meshgrid(lonrange,latrange);
[x,y] = ll2ps(lat,lon);

% extract constituent list
conlist = extract_conlist(Model);
nconst = length(conlist);

% extract top constituent (based on amp) for each grid space
[amp,~,~,~] = tmd_extract_HC(Model,lat,lon,'z');

% reorganize dataset so that each lon and lat corresponds max amp const
dom_const_idx = NaN(length(latrange), length(lonrange));
for i = 1:length(latrange)
    for j = 1:length(lonrange)
         values = amp(:,i,j);

        if all(isnan(values))
            dom_const_idx(i,j) = NaN;
        else
            [~, maxIndex] = max(values, [], 'omitnan');
            dom_const_idx(i,j) = maxIndex;
        end
    end
end

figure
% ice speed plot
ax1 = axes;
measuresps('speed','log')
hold on
measuresps('gl', 'k')

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

% constituent plot
ax2 = axes;

h = pcolor(ax2, x, y, dom_const_idx);
shading(ax2,'flat')

% Make valid constituent cells partially transparent
h.AlphaData = 1 .* ~isnan(dom_const_idx);
h.FaceAlpha = 'flat';
h.AlphaDataMapping = 'none';

cmap = [
    0.1216 0.4667 0.7059
    1.0000 0.4980 0.0549
    0.1725 0.6275 0.1725
    0.8392 0.1529 0.1569
    0.5804 0.4039 0.7412
    0.5490 0.3373 0.2941
    0.8902 0.4667 0.7608
    0.4980 0.4980 0.4980
    0.7373 0.7412 0.1333
    0.0902 0.7451 0.8118
];

colormap(ax2, cmap);
clim(ax2, [0.5, nconst + 0.5]);
cb2 = colorbar(ax2,'eastoutside');
cb2.Ticks = 1:nconst;
cb2.TickLabels = conlist;
cb2.Label.String = 'Dominant tidal constituent';

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
uistack(ax2,'top')

%% create plots showing temporal frequency requirements
% nyquist sampling criteria requires sampling every half period

clear;clc;

% loading CATS2008 Model
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008');
addpath('/Users/jeremywang/Documents/MATLAB/measures_v3.1.2/measures');
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008'; 

% mapping region
latrange = -90.0:0.02:-60;
lonrange = -180:0.02:180;
[lon,lat] = meshgrid(lonrange,latrange);
[x,y] = ll2ps(lat,lon);

% extract constituent list
conlist = extract_conlist(Model);
nconst = length(conlist);

% extract the periods and minimum sampling period
period = zeros(nconst,1);              % hours
freq = zeros(nconst,1);                % cycles/hour
min_sampling_period = zeros(nconst,1); % hours
for i = 1:nconst
    con = conlist{i};
    % extract the period in hours
    [~,~,~,omega,~,~] = constit(con);
    T = 2*pi/omega * (1/(3600)); 
    
    freq(i) = 1/T; % cycles per hours 
    period(i) = T; % hours
    min_sampling_period(i) = T/2; 
end

% extract top constituent (based on amp) for each grid space
[amp,~,~,~] = tmd_extract_HC(Model,lat,lon,'z');

% reorganize dataset so that each lon and lat corresponds max amp const
dom_const_min_period = NaN(length(latrange), length(lonrange));
for i = 1:length(latrange)
    for j = 1:length(lonrange)
         values = amp(:,i,j);

        if all(isnan(values))
            dom_const_min_period(i,j) = NaN;
        else
            [~, maxIndex] = max(values, [], 'omitnan');
            dom_const_min_period(i,j) = min_sampling_period(maxIndex);
        end
    end
end

% plotting values
figure
% ice speed plot
ax1 = axes;
measuresps('speed','log')
hold on
measuresps('gl', 'k')

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

% constituent plot
ax2 = axes;

pcolor(ax2, x, y, dom_const_min_period);
shading(ax2,'flat')

colormap(ax2, winter)
cb2 = colorbar(ax2,'eastoutside');
cb2.Label.String = 'Minimum Sampling Period (hours)';

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
uistack(ax2,'top')

%% Determine the campaign length
% based on the top two amplitudes, calculate the beat frequency

clear;clc;

% loading CATS2008 Model
addpath('/Users/jeremywang/Documents/MATLAB/CATS2008');
addpath('/Users/jeremywang/Documents/MATLAB/measures_v3.1.2/measures');
Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008'; 

% mapping region
latrange = -90.0:0.02:-60;
lonrange = -180:0.02:180;
[lon,lat] = meshgrid(lonrange,latrange);
[x,y] = ll2ps(lat,lon);

% extract constituent list
conlist = extract_conlist(Model);
nconst = length(conlist);

period = struct();
freq = struct();
for i = 1:length(conlist)
    con = conlist{i};
    % extract the period in days
    [~,~,~,omega,~,~] = constit(con);
    T = 2*pi/omega * (1/(3600*24)); % converted to days
    
    freq.(con) = 1/T; % cycles per day 
    period.(con) = T; % days
end

% extract top constituent (based on amp) for each grid space
[amp,~,~,~] = tmd_extract_HC(Model,lat,lon,'z');

% reorganize dataset so that each lon and lat corresponds max amp const
beat_period = NaN(length(latrange), length(lonrange));
beat_freq = NaN(length(latrange), length(lonrange));

for i = 1:length(latrange)
    for j = 1:length(lonrange)
         values = amp(:,i,j);

        if all(isnan(values))
            beat_period(i,j) = NaN;
        else
            [sortedAmp, idx] = sort(values, 'descend');
            con1 = conlist{idx(1)};
            con2 = conlist{idx(2)};
            
            freq1 = freq.(con1);
            freq2 = freq.(con2);
            % calculate the beat period
            beat_freq(i,j) = abs(freq1 - freq2); % cycles/day
            beat_period(i,j) = 1/beat_freq(i,j); % days
        end
    end
end

%%
figure
% ice speed plot
ax1 = axes;
measuresps('speed','log')
hold on
measuresps('gl', 'k')

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

% beat frequency plot
ax2 = axes;

pcolor(ax2, x, y, beat_period);
shading(ax2,'flat')

colormap(ax2, sky)
cb2 = colorbar(ax2,'eastoutside');
cb2.Label.String = 'Beat Period of Top 2 Constituents (Days)';

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
uistack(ax2,'top')
