% =========================================================================
% Table showing decomposition of term premiums variance
% =========================================================================

maturities_in_year = [2;5;10];
H = max(frequency*maturities_in_year);

[uncond_r,uncond_rn,uncond_TPr,uncond_TPrn,uncond_IRP,...
    all_stdv_r,all_stdv_rn,all_stdv_TPr,all_stdv_TPrn,all_stdv_IRP] =...
    compute_uncond_yds_TP(model_sol,H);


model_sol_new = model_sol;
model_sol_new.param_transf.sigma_g

% model.param_transf.mu_gamma0 = RA_values(j);
% model.param   = make_param(model);
% model_sol_new = make_model_sol(model);




