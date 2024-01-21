function[parameters] = make_param(model)
% This function creates a param vector that contains (modifed) parameters.

parameters = zeros(19,1);

parameters(1) = log(model.param_transf.rho_g/(1 - model.param_transf.rho_g));
parameters(2) = log(model.param_transf.rho_z/(1 - model.param_transf.rho_z));
parameters(3) = log(model.param_transf.rho_w/(1 - model.param_transf.rho_w));
parameters(4) = log(model.param_transf.rho_k/(1 - model.param_transf.rho_k));

parameters(5) = log(model.param_transf.sigma_g);
parameters(6) = log(model.param_transf.sigma_z);
parameters(7) = log(model.param_transf.sigma_w);
parameters(8) = log(model.param_transf.stdv_k);

parameters(9) = sign(model.param_transf.mu_gamma1_g) * abs(model.param_transf.mu_gamma1_g)^(1/5);
parameters(10) = sign(model.param_transf.mu_gamma1_z) * abs(model.param_transf.mu_gamma1_z)^(1/5);

parameters(11) = log(model.param_transf.rho_pi/(1 - model.param_transf.rho_pi));

parameters(12) = log(model.param_transf.rho_pi_star/(1 - model.param_transf.rho_pi_star));
parameters(13) = log(model.param_transf.sigma_pi_star);

parameters(14) = log(model.param_transf.sigma_pi_z);

parameters(15) = log(model.param_transf.mu_pi);
parameters(16) = log(model.param_transf.mu_c);
parameters(17) = log(model.param_transf.mu_gamma0);
parameters(18) = log(model.param_transf.delta/(1 - model.param_transf.delta));

parameters(19) = sign(model.param_transf.rho_gz) * 10 * abs(model.param_transf.rho_gz)^(1/3);
