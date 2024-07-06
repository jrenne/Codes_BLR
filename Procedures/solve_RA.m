function[sigma_m] = solve_RA(model)
% This function solves for sigma_m

%indic_w = 4;
indic_m = 5;

%rho_w = model.Phi(indic_w,indic_w);
rho_m = model.Phi(indic_m,indic_m);

%sigma_w = model.Sigma(indic_w,indic_w);
%sigma_m = model.Sigma(indic_m,indic_m); % to be solved for below

% Set unconditional standard deviation for change in slow risk aversion
std_delta_ra = 1;

% Get associated conditional standard deviation of m process
sigma_m = sqrt(1+rho_m)*std_delta_ra/sqrt(2);