function[distance] = compute_distance(sub_parameters,model)
% This functions computes the distance between model-implied moments and targets.

global max_abs_param;

distance = compute_distance_aux(sub_parameters,model);
distance = sum(distance);

% Add penalty for parameter that is too large:
aux = (abs(sub_parameters)>max_abs_param).*...
    (abs(sub_parameters)-max_abs_param);
penalty = 10000 * sum(aux);

distance = distance + penalty;


