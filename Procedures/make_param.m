function[parameters] = make_param(model)
% This function creates a param vector that contains (modified) parameters.

parameters = zeros(20,1);

parameters(1)  = log(model.param_transf.rho_g/(1 - model.param_transf.rho_g));
parameters(2)  = log(model.param_transf.rho_z/(1 - model.param_transf.rho_z));
parameters(3)  = log(model.param_transf.rho_w/(1 - model.param_transf.rho_w));
parameters(4)  = log(model.param_transf.rho_m/(1 - model.param_transf.rho_m));
parameters(5)  = log(model.param_transf.rho_k/(1 - model.param_transf.rho_k));

parameters(6)  = log(model.param_transf.sigma_g);
parameters(7)  = log(model.param_transf.sigma_z);
parameters(8)  = log(model.param_transf.sigma_w);
parameters(9)  = log(model.param_transf.sigma_m);
parameters(10) = log(model.param_transf.sigma_k);

parameters(11) = log(model.param_transf.rho_pi_star/(1 - model.param_transf.rho_pi_star));
parameters(12) = log(model.param_transf.sigma_pi_star);

parameters(13) = log(model.param_transf.rho_pi_tilde/(1 - model.param_transf.rho_pi_tilde));
parameters(14) = log(model.param_transf.sigma_pi_z);

parameters(15) = log(model.param_transf.mu_pi);
parameters(16) = log(model.param_transf.mu_c);
parameters(17) = log(model.param_transf.mu_gamma);
parameters(18) = log((1+model.param_transf.mu_kappa)/(1 - model.param_transf.mu_kappa));
parameters(19) = log(model.param_transf.delta/(1 - model.param_transf.delta));

parameters(20) = log(model.param_transf.rho_gz);