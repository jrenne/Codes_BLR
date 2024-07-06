function[model_sol] = make_model_sol(model,endo_flag)
% This functions defines the model.
% endo_flag: calculate mu_c/mu_gamma endogenously (default: 1)

model_sol = model;

parameters = model.param;

model_sol.names_param = {'rho_g';'rho_z';'rho_w';'rho_m';'rho_k';...
    'sigma_g';'sigma_z';'sigma_w';'sigma_m';'sigma_k';...
    'rho_pi_star';'sigma_pi_star';'rho_pi_tilde';'sigma_pi_z';'mu_pi';'mu_c';...
    'mu_gamma';'mu_kappa';'delta';'rho_gz'};
model_sol.names_param_Latex = {'$\rho_g$';'$\rho_z$';'$\rho_w$';'$\rho_m$';'$\rho_k$';...
    '$\sigma_g$';'$\sigma_z$';'$\sigma_w$';'$\sigma_m$';'$\sigma_k$';...
    '$\rho_{\pi^*}$';'$\sigma_{\pi^*}$';'$\rho_{\tilde{\pi}}$';'$\sigma_{\pi,z}$';'$\mu_\pi$';'$\mu_c$';...
    '$\mu_\gamma$';'$\mu_\kappa$';'$\delta$';'$\rho_{g,z}$'};
model_sol.names_shocks_Latex = {'$\varepsilon_{g,t}$';'$\varepsilon_{z,1,t}$';'$\varepsilon_{z,2,t}$';...
    '$\varepsilon_{w,t}$';'$\varepsilon_{m,t}$';'$\varepsilon_{k,t}$';'$\sigma_{\pi^*,t}$'};

[rho_g, rho_z, rho_w, rho_m, rho_k,...
    sigma_g, sigma_z, sigma_w, sigma_m, sigma_k,...
    rho_pi_star, sigma_pi_star, rho_pi_tilde, sigma_pi_z, mu_pi, mu_c,...
    mu_gamma, mu_kappa, delta, rho_gz] = make_param_transf(parameters);

% X_t is as follows:
%     [g_t,z_t,t_{t-1},w_t,m_t,k_t]'
% eps_t is as follows:
%     [eps_{g,t},eps_{z,1,t},eps_{z,2,t},eps_{w,t},eps_{m,t},eps_{k,t},eps_{pi*,t}]'

n_X   = 6;
n_Z   = 2;
n_eps = 7;
n_Y   = n_X + n_Z + n_X^2;

model_sol.mu_c0     = mu_c;
model_sol.mu_c1     = zeros(n_X,1);
model_sol.mu_c1(1)  = 1;
model_sol.mu_c1(2)  = 1;
model_sol.mu_c1(3)  = -1;

model_sol.delta     = delta;

model_sol.mu_gamma0    = mu_gamma;
model_sol.mu_gamma1    = zeros(n_X,1);
model_sol.mu_gamma1(4) = 1; % w
model_sol.mu_gamma1(5) = 1; % m

model_sol.Phi       = 0 * eye(n_X);
model_sol.Phi(1,1)  = rho_g;
model_sol.Phi(2,2)  = rho_z;
model_sol.Phi(3,2)  = 1;
model_sol.Phi(4,4)  = rho_w;
model_sol.Phi(5,5)  = rho_m;
model_sol.Phi(6,6)  = rho_k;
model_sol.Phi(1,2)  = rho_gz;

model_sol.Sigma      = zeros(n_X,n_eps);
model_sol.Sigma(1,1) = sigma_g;
model_sol.Sigma(2,2) = 1/sqrt(2)*sigma_z;
model_sol.Sigma(2,3) = 1/sqrt(2)*sigma_z;
model_sol.Sigma(4,4) = sigma_w;
model_sol.Sigma(5,5) = sigma_m;
model_sol.Sigma(6,6) = sigma_k;

model_sol.Phi_Z      = zeros(n_Z,n_Z);
model_sol.Phi_Z(1,1) = rho_pi_star;
model_sol.Phi_Z(2,2) = rho_pi_tilde;

model_sol.mu_kappa   = mu_kappa;

model_sol.Gamma0     = zeros(n_Z*n_eps,1);
model_sol.Gamma0(4)  = +(1+mu_kappa)*sigma_pi_z/2;
model_sol.Gamma0(6)  = -(1-mu_kappa)*sigma_pi_z/2;
model_sol.Gamma0(13) = sigma_pi_star;

model_sol.Gamma1      = zeros(n_Z*n_eps,n_X);
model_sol.Gamma1(4,6) = sigma_pi_z/2;
model_sol.Gamma1(6,6) = sigma_pi_z/2;

% Define inflation process:
model_sol.mu_pi0    = mu_pi;
model_sol.mu_piZ    = [1;1];
model_sol.mu_piX    = zeros(n_X,1);
model_sol.mu_piXX   = zeros(n_X*(n_X+1)/2,1);

if nargin == 1
    endo_flag = 1;
elseif ~ismember(endo_flag,[0 1])
    disp("Invalid value for endo_flag (default of '1' is used).");
    endo_flag = 1;
end

if endo_flag == 1

    % Calculate endogenous parameters:
    % Ordering such that first those endogenous parameters are calculated
    % which are independent of the other endogenous parameters.

    % overwrite sigma_m to constrain risk aversion process
    sigma_m = solve_RA(model_sol);
    indic_param = strcmp(model_sol.names_param,'sigma_m');
    model_sol.param(indic_param) = log(sigma_m);
    indic_m = 5;
    model_sol.Sigma(indic_m,indic_m) = sigma_m;

    % overwrite sigma_k to constrain kappa process
    sigma_k = solve_kappa(model_sol);
    indic_param = strcmp(model_sol.names_param,'sigma_k');
    model_sol.param(indic_param) = log(sigma_k);
    indic_k = 6;
    model_sol.Sigma(indic_k,indic_k) = sigma_k;

    % overwrite mu_c0 and mu_gamma0 to match real rates
    [mu_c,mu_gamma]     = solve_rr(model_sol);
    indic_param = strcmp(model_sol.names_param,'mu_c');
    model_sol.param(indic_param) = log(mu_c);
    indic_param = strcmp(model_sol.names_param,'mu_gamma');
    model_sol.param(indic_param) = log(mu_gamma);
    model_sol.mu_c0     = mu_c;
    model_sol.mu_gamma0 = mu_gamma;

end

% Solve model:
model_sol = compute_sdf(model_sol);

model_sol.param_transf = [rho_g rho_z rho_w rho_m rho_k...
    sigma_g sigma_z sigma_w sigma_m sigma_k...
    rho_pi_star sigma_pi_star rho_pi_tilde sigma_pi_z mu_pi mu_c...
    mu_gamma mu_kappa delta rho_gz];
