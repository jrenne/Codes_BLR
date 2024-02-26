function[model_sol] = make_model_sol(model)
% This functions computes the distance between model-implied moments and targets.

model_sol = model;

parameters = model.param;

rho_g   = .999 * exp(parameters(1))/(1+exp(parameters(1)));
rho_z   = exp(parameters(2))/(1+exp(parameters(2)));
rho_w   = exp(parameters(3))/(1+exp(parameters(3)));
rho_k   = exp(parameters(4))/(1+exp(parameters(4)));

sigma_g = exp(parameters(5));

% =========================================================================
sigma_z = exp(parameters(6)) * sqrt(1-rho_z^2);
% =========================================================================

sigma_w = exp(parameters(7));
stdv_k  = exp(parameters(8));
sigma_k = stdv_k * sqrt(1 - rho_k^2);

mu_gamma1_g = - exp(parameters(9));
mu_gamma1_z = - exp(parameters(10));

rho_pi   = exp(parameters(11))/(1+exp(parameters(11)));

rho_pi_star   = exp(parameters(12))/(1+exp(parameters(12)));
sigma_pi_star = exp(parameters(13));

sigma_pi_z = exp(parameters(14));

mu_pi     = exp(parameters(15));
mu_c      = exp(parameters(16));
mu_gamma0 = exp(parameters(17));
delta     = exp(parameters(18))/(1+exp(parameters(18)));

rho_gz = exp(parameters(19));

model_sol.param_transf = [rho_g rho_z rho_w rho_k...
    sigma_g sigma_z sigma_w sigma_k...
    mu_gamma1_g mu_gamma1_z rho_pi rho_pi_star...
    sigma_pi_star sigma_pi_z mu_pi mu_c...
    mu_gamma0 delta rho_gz];
model_sol.names_param = {'rho_g';'rho_z';'rho_w';'rho_k';...
    'sigma_g';'sigma_z';'sigma_w';'sigma_k';...
    'mu_gamma1_g';'mu_gamma1_z';'rho_pi';'rho_pi_star';...
    'sigma_pi_star';'sigma_pi_z';'mu_pi';'mu_c';...
    'mu_gamma0';'delta';'rho_gz'};
model_sol.names_param_Latex = {'$\rho_g$';'$\rho_z$';'$\rho_w$';'$\rho_k$';...
    '$\sigma_g$';'$\sigma_z$';'$\sigma_w$';'$\sigma_k$';...
    '$\mu_{\gamma,1,g}$';'$\mu_{\gamma,1,z}$';'$\rho_\pi$';'$\rho_{\pi}^*$';...
    '$\sigma_{\pi}^*$';'$\sigma_{\pi,z}$';'$\mu_\pi$';'$\mu_c$';...
    '$\mu_{\gamma,0}$';'$\delta$';'$\rho_{g,z}$'};

% X_t is as follows:
%     [g_t,z_t,t_{t-1},w_t,k_t]'
% eps_t is as follows:
%     [eps_{g,t},eps_{z,1,t},eps_{z,2,t},eps_{w,t},eps_{k,t},eps_{pi*,t}]'

model_sol.names_shocks_Latex = {'$\varepsilon_{g,t}$';'$\varepsilon_{z,1,t}$';...
    '$\varepsilon_{z,2,t}$';...
    '$\varepsilon_{w,t}$';'$\varepsilon_{k,t}$';'$\sigma_{\pi^*,t}$'};

n_X   = 5;
n_Z   = 2;
n_eps = 6;
n_Y   = n_X + n_Z + n_X^2;

model_sol.mu_c0     = mu_c;
model_sol.mu_c1     = zeros(n_X,1);
model_sol.mu_c1(1)  = 1;
model_sol.mu_c1(2)  = 1;
model_sol.mu_c1(3)  = -1;

model_sol.delta     = delta;

model_sol.mu_gamma0 = mu_gamma0;

model_sol.mu_gamma1    = zeros(n_X,1);
model_sol.mu_gamma1(1) = mu_gamma1_g;
model_sol.mu_gamma1(2) = mu_gamma1_z;
model_sol.mu_gamma1(4) = 1; % w

model_sol.Phi       = 0 * eye(n_X);
model_sol.Phi(1,1)  = rho_g;
model_sol.Phi(2,2)  = rho_z;
model_sol.Phi(3,2)  = 1;
model_sol.Phi(4,4)  = rho_w;
model_sol.Phi(5,5)  = rho_k;
model_sol.Phi(1,2)  = rho_gz;

model_sol.Sigma      = zeros(n_X,n_eps);
model_sol.Sigma(1,1) = sigma_g;
model_sol.Sigma(2,2) = 1/sqrt(2)*sigma_z;
model_sol.Sigma(2,3) = 1/sqrt(2)*sigma_z;
model_sol.Sigma(4,4) = sigma_w;
model_sol.Sigma(5,5) = sigma_k;

model_sol.Phi_Z      = zeros(2,2);
model_sol.Phi_Z(1,1) = rho_pi_star;
model_sol.Phi_Z(2,2) = rho_pi;

model_sol.Gamma0     = zeros(n_Z*n_eps,1);
model_sol.Gamma0(4)  = +sigma_pi_z/2;
model_sol.Gamma0(6)  = -sigma_pi_z/2;
model_sol.Gamma0(11) = sigma_pi_star;

model_sol.Gamma1      = zeros(n_Z*n_eps,n_X);
model_sol.Gamma1(4,5) = sigma_pi_z/2;
model_sol.Gamma1(6,5) = sigma_pi_z/2;

% Define inflation process:
model_sol.mu_pi0    = mu_pi;
model_sol.mu_piZ    = [1;1];
model_sol.mu_piX    = zeros(n_X,1);
model_sol.mu_piXX   = zeros(n_X*(n_X+1)/2,1);

% Solve model:
model_sol = compute_sdf(model_sol);



