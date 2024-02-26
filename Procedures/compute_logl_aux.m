function[logl,A4r,B4r,A4rn,B4rn,C4rn,D4rn,model_sol_new] = compute_logl_aux(sub_parameters,Data_StateSpace,model)
% This functions computes the distance between model-implied moments and targets.

global FILTER;

new_parameters = model.param;
new_parameters(FILTER==1) = sub_parameters;

model_new = model;
model_new.param = new_parameters;

model_sol_new = make_model_sol(model_new);

[StateSpace,xi_00,P_00,A4r,B4r,A4rn,B4rn,C4rn,D4rn] = prepare_State_Space(model_sol_new,Data_StateSpace);

n_X = size(model_sol_new.PhiQ,1);

[~,~,logl] = kalman(StateSpace,Data_StateSpace.dataset,xi_00,P_00,n_X);



