function[distance] = compute_logl(sub_parameters,Data_StateSpace,model,indic_add_moments)
% This functions computes the distance between model-implied moments and targets.

global max_abs_param moments;

[logl,A4r,B4r,A4rn,B4rn,C4rn,D4rn,model_sol_new] = compute_logl_aux(sub_parameters,Data_StateSpace,model);

% Add moment-based term in distance
if indic_add_moments == 1

    mim = compute_moments(model_sol_new,A4r,B4r,A4rn,B4rn,C4rn,D4rn);
    targets = moments.values;
    distance = mim - targets;
    
    % Special treatment for min and max:
    indic_set2zero = find((distance>0).*(moments.indic_min==1));
    distance(indic_set2zero) = 0;
    indic_set2zero = find((distance<0).*(moments.indic_max==1));
    distance(indic_set2zero) = 0;

    distance = 100*(distance.^2).*moments.weights;
    distance = sum(distance);
else
    distance = 0;
end

% Add penalty for parameter that is too large:
aux = (abs(sub_parameters)>max_abs_param).*...
    (abs(sub_parameters)-max_abs_param);
penalty = 10000 * sum(aux);

distance = - logl +  100*distance + penalty;

