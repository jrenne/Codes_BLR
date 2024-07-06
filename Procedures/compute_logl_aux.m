function[logl,logl_t,A4r,B4r,A4rn,B4rn,C4rn,D4rn,model_sol_new] = compute_logl_aux(sub_parameters,Data_StateSpace,model,endo_flag)
% This functions computes the log-likelihood

global FILTER;

new_parameters = model.param;
new_parameters(FILTER==1) = sub_parameters;

model_new = model;
model_new.param = new_parameters;

if nargin == 3
    endo_flag = 1;
elseif ~ismember(endo_flag,[0 1])
    disp("Invalid value for endo_flag (default of '1' is used).");
    endo_flag = 1;
end

model_sol_new = make_model_sol(model_new,endo_flag);

[StateSpace,xi_00,P_00,A4r,B4r,A4rn,B4rn,C4rn,D4rn] = prepare_State_Space(model_sol_new,Data_StateSpace);

n_X = size(model_sol_new.PhiQ,1);

[~,~,logl,logl_t] = kalman(StateSpace,Data_StateSpace.dataset,xi_00,P_00,n_X);