% ========== Define initial model =========================================

model = struct;
model.param_transf = struct;

model.param_transf.mu_c    = mean(data_bis(:,find(strcmp(data_names,{'real GDP'}))),'omitnan')/100;

model.param_transf.rho_g   = (.995)^(1/frequency); %.996118;
model.param_transf.rho_gz  = 0.01/frequency; %.0112;
model.param_transf.sigma_g = 0.001*sqrt(1-model.param_transf.rho_g^2); %0.000001*sqrt(1-rho_g^2);

model.param_transf.rho_d   = (.6)^(1/frequency); %.3815;
model.param_transf.sigma_d = 0.02/frequency*sqrt(1-model.param_transf.rho_d^2); %0.0207*sqrt(1-rho_z^2);

model.param_transf.rho_s   = (.6)^(1/frequency); %.3815;
model.param_transf.sigma_s = 0.02/frequency*sqrt(1-model.param_transf.rho_s^2); %0.0207*sqrt(1-rho_z^2);

model.param_transf.rho_w   = .98;
model.param_transf.sigma_w = 3 * sqrt(1-model.param_transf.rho_w^2); %.1*sqrt(1-rho_w^2);

model.param_transf.rho_k   = (.98)^(1/frequency);
model.param_transf.stdv_k  = .3;

model.param_transf.delta = (.995)^(1/frequency);

model.param_transf.mu_gamma0   = 10;
model.param_transf.mu_gamma1_g = 0.00001;
model.param_transf.mu_gamma1_z = -100;

model.param_transf.mu_pi = mean(data_bis(:,find(strcmp(data_names,{'CPI all '}))),'omitnan')/100;

model.param_transf.rho_pi   = .9;
model.param_transf.sigma_pi = 0.0001*sqrt(1-model.param_transf.rho_pi^2);

model.param_transf.rho_pi_star   = .98;
model.param_transf.sigma_pi_star = 0.01*sqrt(1-model.param_transf.rho_pi_star^2);

model.param_transf.sigma_pi_g = 0.01*sqrt(1-model.param_transf.rho_pi^2);
model.param_transf.sigma_pi_d = 0.002*sqrt(1-model.param_transf.rho_pi^2);
model.param_transf.sigma_pi_s = 0.002*sqrt(1-model.param_transf.rho_pi^2);

model.param = make_param(model);



