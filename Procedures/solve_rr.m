function[mu_c0_sol,mu_gamma0_sol] = solve_rr(model)
% This function solves for average consumption growth and risk aversion
% based on the other parameters of the model to match the real rate curve.

global frequency;

% Load model

n_X = size(model.Phi,1);

%mu_c0     = model.mu_c0;     % to be solved for below
mu_c1     = model.mu_c1;
delta     = model.delta;
%mu_gamma0 = model.mu_gamma0; % to be solved for below
mu_gamma1 = model.mu_gamma1;
Phi       = model.Phi;
Sigma     = model.Sigma;

% declare symbolic math parameters

syms 'mu_c0' 'mu_gamma0' real;

% calculate real SDF as a function of the symbolic math parameters

lambda0 = ((1-mu_gamma0)*delta*(eye(n_X)-delta*Phi')^(-1)*Phi' - mu_gamma0*eye(n_X)) * mu_c1;
lambda1 = -(mu_gamma1*mu_c1')*(delta*Phi*(eye(n_X)-delta*Phi)^(-1)+eye(n_X));
eta0 = -log(delta) + mu_c0 + 0.5*mu_c1'*Sigma*(Sigma')*mu_c1 + lambda0'*Sigma*(Sigma')*mu_c1;
eta1 = (lambda1*Sigma*(Sigma')+Phi')*mu_c1;
muQ = Sigma*(Sigma')*lambda0;
PhiQ = Phi + Sigma*(Sigma')*lambda1';

% Calculate unconditional moments of real rates as a function of the symbolic math parameters

ah = -eta0;
bh = -eta1;

H = frequency * 10;
rr = nan(H,1);
rr = sym(rr);
rr(1) = 100 * frequency * eta0;
for h = 2:H
    ah = -eta0 + ah + bh'*muQ + 0.5*bh'*Sigma*(Sigma')*bh;
    bh = -eta1 + PhiQ'*bh;
    rr(h) = -ah/h * 100 * frequency;
end

% Solve for mu_c0 & mu_gamma0 to match calibrated real rates

mat_short_years = 2;
mat_long_years  = 10;

rr_cal(1) = 1.5;
rr_cal(2) = 2.3;
[mu_c0_sol,mu_gamma0_sol] = solve(rr(mat_short_years*frequency) == rr_cal(1),rr(mat_long_years*frequency) == rr_cal(2));

mu_c0_sol     = double(mu_c0_sol);
mu_gamma0_sol = double(mu_gamma0_sol);