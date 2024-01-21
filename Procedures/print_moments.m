function [M] = print_moments(sub_parameters,model_sol,indic_disp)
% print model-implied moments

global FILTER moments frequency;

new_parameters = model_sol.param;
new_parameters(FILTER==1) = sub_parameters;

model_new = model_sol;
model_new.param = new_parameters;
model_sol_new = make_model_sol(model_new);

% Compute factor loadings:
H = 10  * frequency; % H is expressed in quarters
[~,~,A4r,B4r] = compute_AB(model_sol,H);
[~,~,~,~,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol_new,H);

% Compare observed and model-implied moments (at initial parameters):
M = array2table(round([moments.values compute_moments(model_sol_new,A4r,B4r,A4rn,B4rn,C4rn,D4rn)...
    compute_distance_aux(sub_parameters,model_sol_new)...
    moments.weights moments.indic_min moments.indic_max],4),...
    'RowNames',moments.names,...
    'VariableNames',{'Data','Model','Loss contrib.','Weight','Min','Max'});

if(indic_disp==1)
    disp(M);
end
