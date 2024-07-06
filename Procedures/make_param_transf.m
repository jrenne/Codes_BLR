function[rho_g, rho_z, rho_w, rho_m, rho_k,...
    sigma_g, sigma_z, sigma_w, sigma_m, sigma_k,...
    rho_pi_star, sigma_pi_star, rho_pi_tilde, sigma_pi_z, mu_pi, mu_c,...
    mu_gamma, mu_kappa, delta, rho_gz] = make_param_transf(parameters)

rho_g   = exp(parameters(1))/(1+exp(parameters(1)));
rho_z   = exp(parameters(2))/(1+exp(parameters(2)));
rho_w   = exp(parameters(3))/(1+exp(parameters(3)));
rho_m   = exp(parameters(4))/(1+exp(parameters(4)));
rho_k   = exp(parameters(5))/(1+exp(parameters(5)));

sigma_g = exp(parameters(6));
sigma_z = exp(parameters(7));
sigma_w = exp(parameters(8));
sigma_m = exp(parameters(9));
sigma_k = exp(parameters(10));

rho_pi_star   = exp(parameters(11))/(1+exp(parameters(11)));
sigma_pi_star = exp(parameters(12));

rho_pi_tilde = exp(parameters(13))/(1+exp(parameters(13)));
sigma_pi_z   = exp(parameters(14));

mu_pi    = exp(parameters(15));
mu_c     = exp(parameters(16));
mu_gamma = exp(parameters(17));
mu_kappa = (exp(parameters(18))-1)/(exp(parameters(18))+1);
delta    = exp(parameters(19))/(1+exp(parameters(19)));

rho_gz = exp(parameters(20));
