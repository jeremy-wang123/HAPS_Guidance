function amp_cons = extract_mean_amp(region)
    % extracting the mean amplitude of constituents in a region
    addpath('/Users/jeremywang/Documents/MATLAB/CATS2008')
    Model = '/Users/jeremywang/Documents/MATLAB/CATS2008/Model_CATS2008';
    
    fig = figure('Visible','off');
    
    % extracting region of interest
    mapzoomps(region)
    xl = xlim;
    yl = ylim;
    
    x_corner = [xl(1) xl(2) xl(2) xl(1)];
    y_corner = [yl(1) yl(1) yl(2) yl(2)];
    
    [lat_corner, lon_corner] = ps2ll(x_corner, y_corner);
    
    latlim = [min(lat_corner) max(lat_corner)];
    lonlim = [min(lon_corner) max(lon_corner)];
    
    latrange = latlim(1):0.02:latlim(2);
    lonrange = lonlim(1):0.02:lonlim(2);
    
    % Get elevation file from model control file
    [ModName, GridName, Fxy_ll] = rdModFile(Model, 1);
    
    conList = rd_con(ModName);
    conList_clean = lower(strtrim(cellstr(conList)));
    
    [lat,lon] = meshgrid(latrange,lonrange);
    
    [amps,~,~,~] = tmd_extract_HC(Model,lat,lon,'z');
    means = squeeze(mean(amps, [2 3], 'omitnan'));
    
    cons = extract_conlist(Model);
    amp_cons = struct();
    for i=1:length(cons)
       con = cons{i};
       amp_cons.(con) = means(i);
    end
end