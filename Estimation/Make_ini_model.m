% ========== Define initial model =========================================

model = struct;
model.param_transf = struct;

model.param_transf.rho_g   = 0.99;
model.param_transf.sigma_g = 0.0005;

model.param_transf.rho_z   = 0.95;
model.param_transf.sigma_z = 0.05;

model.param_transf.rho_w   = 0.9;
model.param_transf.sigma_w = 5;

model.param_transf.rho_m   = 0.99;
model.param_transf.sigma_m = 1;

model.param_transf.rho_k   = 0.95;
model.param_transf.sigma_k = 0.1;

model.param_transf.delta   = 0.98;
model.param_transf.rho_gz  = 0.0005;

model.param_transf.mu_gamma   = 20;

model.param_transf.mu_kappa = -0.3;

model.param_transf.mu_pi = mean(data(:,find(strcmp(data_names,{'CPI inflation'}))),'omitnan')/100;
model.param_transf.mu_c  = mean(data(:,find(strcmp(data_names,{'Consumption growth'}))),'omitnan')/100;

model.param_transf.rho_pi_star   = 0.99;
model.param_transf.sigma_pi_star = 0.0005;

model.param_transf.rho_pi_tilde   = 0.9;
model.param_transf.sigma_pi_z = 0.01;

model.param = make_param(model);
