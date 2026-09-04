function near_grounding_line = is_near_grounding_line(x0, y0, radius, spacing)
%IS_NEAR_GROUNDING_LINE Checks whether grounding status changes nearby.
%
% near_grounding_line = is_near_grounding_line(x0, y0, radius, spacing)
%
% Inputs:
%   x0, y0  - Center coordinates in meters
%   radius  - Search radius in meters, e.g. 10e3
%   spacing - Sampling-grid spacing in meters, e.g. 500
%
% Output:
%   near_grounding_line - true if any sampled point within the radius has
%                         a different grounding status than the center

    % Grounding status at the center
    center_status = isgrounded(x0, y0);

    % Create a local square grid around the center
    offsets = -radius:spacing:radius;
    [dx, dy] = meshgrid(offsets, offsets);

    % Keep only points inside the circular search region
    inside_circle = hypot(dx, dy) <= radius;

    x_test = x0 + dx(inside_circle);
    y_test = y0 + dy(inside_circle);

    % Evaluate grounding status at all nearby points
    nearby_status = isgrounded(x_test, y_test);

    % True if at least one nearby point differs from the center
    near_grounding_line = any(nearby_status ~= center_status);
end