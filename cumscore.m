function score = cumscore(distance, angle_spacing, x0, y0)
%CUMSCORE averages the top 20% of weight scores with a 360 degree radius
angles = 0:angle_spacing:179;

weights = zeros(1, length(angles));
for i=1:length(angles)
    weights(i) = extract_weights(distance,angles(i),x0,y0);
end
top = maxk(weights, round(0.2*numel(weights)));
score = mean(top);
end