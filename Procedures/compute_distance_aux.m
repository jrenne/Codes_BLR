function[distance] = compute_distance_aux(sub_parameters,model)
% This functions computes the distance between model-implied moments and targets.

global FILTER moments frequency;

new_parameters = model.param;
new_parameters(FILTER==1) = sub_parameters;

model_new = model;
model_new.param = new_parameters;

model_sol_new = make_model_sol(model_new);

% Compute factor loadings:
H = 10  * frequency; % H is expressed in quarters
[~,~,A4r,B4r] = compute_AB(model_sol_new,H);
[~,~,~,~,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol_new,H);

% Compute moments:
mim = compute_moments(model_sol_new,A4r,B4r,A4rn,B4rn,C4rn,D4rn);

targets = moments.values;

distance = mim - targets;

% Special treatment for min and max:
indic_set2zero = find((distance>0).*(moments.indic_min==1));
distance(indic_set2zero) = 0;
indic_set2zero = find((distance<0).*(moments.indic_max==1));
distance(indic_set2zero) = 0;

distance = 100*(distance.^2).*moments.weights;


