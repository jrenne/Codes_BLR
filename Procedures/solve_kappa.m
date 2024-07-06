function[sigma_k] = solve_kappa(model)
% This function solves for sigma_k

indic_k = 6;

mu_kappa = model.mu_kappa;
rho_k      = model.Phi(indic_k,indic_k);
%sigma_k    = model.Sigma(indic_k,indic_k); % to be solved for below

% Set probability for abs(kappa) > bound
alpha = 1/1000;
bound = 1;

% Get associated unconditional standard deviation
if mu_kappa > 0     %upper bound more likely
    z_score = norminv(1-alpha,0,1);
    std_k = (+bound - mu_kappa)/z_score;
else                %lower bound more likely
    z_score = norminv(alpha,0,1);
    std_k = (-bound - mu_kappa)/z_score;
end

% Get associated conditional standard deviation
sigma_k = std_k*sqrt(1-rho_k^2);