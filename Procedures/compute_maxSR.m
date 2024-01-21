function [maxSR] = compute_maxSR(model_sol,Y)
% compute maximum Sharpe ratio

T = size(Y,1);
n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

lambda0 = model_sol.lambda0;
lambda1 = model_sol.lambda1;
maxSR   = zeros(T,1);

for t = 1:T
    u1 = zeros(n_Y,1);
    u2 = zeros(n_Y,1);
    u1(1:n_X) = lambda0    + lambda1 * Y(t,1:n_X)';
    u2(1:n_X) = 2*(lambda0 + lambda1 * Y(t,1:n_X)');

    Psi1 = compute_LT_Y(model_sol,u1);
    Psi2 = compute_LT_Y(model_sol,u2);
    
    psi1 = Psi1.Psi0 + Psi1.Psi_Y' * Y(t,:)';
    psi2 = Psi2.Psi0 + Psi2.Psi_Y' * Y(t,:)';

    maxSR(t) = sqrt(exp(psi2 - 2 * psi1) - 1);
end

