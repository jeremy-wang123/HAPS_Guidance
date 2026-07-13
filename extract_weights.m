function total_weight = extract_weights(distance, theta, x0, y0)
%EXTRACT_WEIGHTS extracts the weight for a given flight
% Inputs: distance (how far the plane flies), theta (angle of flight in
% degrees from the positive x direction), x0 and y0 (starting points in
% polar stereographic coords)
% Output: weighted score of that trajectory

% determine the endpoints of the flight
x1 = x0 + distance*cosd(theta);
y1 = y0 + distance*sind(theta);

% define the line trajectory
x_traj =[x0 x1];
y_traj = [y0 y1];

% extracting the flow lines from the thwaites region
[xstart, ystart] = psgrid('thwaites glacier', 500, 30, 'xy');
[xf,yf] = flowline(xstart, ystart);

%%% Extract the intersections between a given flight path and flow lines
intersections = struct([]);
count = 1;

for f = 1:numel(xf)

    % Extract one flowline
    lat = xf{f};
    lon = yf{f};
    
    % eliminating flow lines smaller than 2 data points
    valid = isfinite(lat) & isfinite(lon);
    lat = lat(valid);
    lon = lon(valid);

    if numel(lat) < 2
        continue
    end

    % Convert flowline lat/lon to x/y meters
    [xline, yline] = ll2ps(lat, lon);

    % Find intersections between candidate line and this flowline
    [xi, yi] = polyxpoly(x_traj, y_traj, xline, yline);

    if isempty(xi)
        continue
    end

    % Store every intersection
    for q = 1:numel(xi)

        % intersections(count).path_id = p; % storing the specific line
        intersections(count).flowline_id = f; % storing the flowline index
        % storing intersection coordinates
        intersections(count).x_intersect = xi(q);
        intersections(count).y_intersect = yi(q);

        % nearest sampled flowline point
        d = hypot(xline - xi(q), yline - yi(q)); % computes distance to every sampled point
        [~, idx] = min(d); % returns index for minimum point

        intersections(count).flowline_index = idx; % storing that index

        count = count + 1;
    end
end

%%% calculate the corresponding weights for each intersection, and sum
% Weighting scheme: 
% since we want to use flowline intersection points as the metric for
% quality of the line, we can segment each individual flow line based on
% their distance to the grounding like on a scale of 0 to 1, with an
% incrementing decrease in weight until we reach 50km. This decrease can be
% linear for now, representing the value of the information collected at a
% given location.

total_weight = 0; 

for i=1:length(intersections)
    % extract key parameters from struct
    flowline_id = intersections(i).flowline_id;
    flowline_index = intersections(i).flowline_index;
    
    % extract lat and lon vectors for the flowline
    lat = xf{flowline_id};
    lon = yf{flowline_id};

    valid = isfinite(lat) & isfinite(lon);
    lat = lat(valid);
    lon = lon(valid);
    
    % grounding mask
    grounded = isgrounded(lat, lon);
    
    % find index of first floating data point
    firstFloating = find(~grounded, 1, 'first');
    
    % skip any flow line that does not reach the grounding line
    if isempty(firstFloating)
        continue;
    end
    
    % index of the grounding line
    gl_index = firstFloating - 1;

    % skip intersections on floating ice
    if flowline_index > gl_index
        intersections(i).gl_distance = NaN;
        intersections(i).weight = 0;
        continue
    end
    
    
    % convert lat and lon into polar
    [x,y] = ll2ps(lat, lon);
    
    % distance between consecutive points along the flowline
    ds = hypot(diff(x), diff(y));
    
    % cumulative distance from beginning of flowline
    cumdist = [0; cumsum(ds)];
    
    % distance along flowline between intersection point and grounding line (km)
    dist_intersection_to_GL = abs(cumdist(gl_index) - cumdist(flowline_index))/1000;
    
    intersections(i).gl_distance = dist_intersection_to_GL; % km

    % linear decrease in weighting beyond 10km to 50km for weights from 0 to 1
    if dist_intersection_to_GL <= 10
        intersections(i).weight = 1;
    elseif dist_intersection_to_GL  >= 50
        intersections(i).weight = 0;
    else
        intersections(i).weight = 1 - (dist_intersection_to_GL - 10)/(50 - 10);
    end

    total_weight = total_weight + intersections(i).weight;
end

end