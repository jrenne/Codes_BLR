% ========== Define initial model =========================================

model = struct;
model.param_transf = struct;


model.param_transf.rho_g   = 0.9998; %.996118;
model.param_transf.rho_gz  = 0.0005; %.0112;
model.param_transf.sigma_g = 0.0000001;

model.param_transf.rho_z   = .9503; %.3815;
model.param_transf.sigma_z = sqrt(2)*.0053; %0.0207*sqrt(1-rho_z^2);

model.param_transf.rho_w   = 0.9889;
model.param_transf.sigma_w = 0.9642; %.1*sqrt(1-rho_w^2);

model.param_transf.rho_k   = 0.9405;
model.param_transf.stdv_k  = .3;

model.param_transf.delta = 0.9989;

model.param_transf.mu_gamma0   = 19.42;
model.param_transf.mu_gamma1_g = 0.000001;
model.param_transf.mu_gamma1_z = 0.000001;

model.param_transf.mu_pi = mean(data_bis(:,find(strcmp(data_names, ...
    {'CPI all '}))),'omitnan')/100;
model.param_transf.mu_c  = mean(data_bis(:,find(strcmp(data_names, ...
    {'Real consumption of all goods and services'}))),'omitnan')/100;

model.param_transf.rho_pi   = 0.9725;

model.param_transf.rho_pi_star   = 0.8678;
model.param_transf.sigma_pi_star = 0.0029;

model.param_transf.sigma_pi_z = 0.002;

model.param = make_param(model);



