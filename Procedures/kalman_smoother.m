function[all_xi_tT,all_P_tT] = kalman_smoother(StateSpace,DataObs,xi_00,P_00,option)

% =========================================================================
% State-space model:
%
% Y(t)  = A' + H'*xi(t) + eps(t)
% xi_t  = mu_xi + F*xi(t-1) + nu(t)
%
% with:
%
% Var(eps(t)) = R
% Var(nu(t)|xi(t-1)) = functionQ(xi(t-1))
%
% We also have vec[Var(nu(t)|xi(t-1))] = Theta0 + Theta1 xi(t-1)
% (this is what functionQ computes)
% =========================================================================


[all_xi_tt,all_P_tt] = kalman(StateSpace,DataObs,xi_00,P_00);

mu_xi = StateSpace.mu_xi;
F     = StateSpace.F;
n_xi  = size(F,1);
T     = size(DataObs,1); % number of periods

% Initialize matrices:
all_xi_tT = all_xi_tt;
all_P_tT  = all_P_tt;

xi_tT = all_xi_tt(T,:)';
P_tT  = reshape(all_P_tT(T,:),n_xi,n_xi);

for T_t = 1:(T-1)

    t = T - T_t;

    xi_tt = all_xi_tt(t,:)';
    P_tt = reshape(all_P_tt(t,:),n_xi,n_xi);

    Q = functionQ(StateSpace,xi_tt);
    P_tp1_t = F * P_tt * F' + Q;
    P_star_t = P_tt * F' * P_tp1_t^(-1);

    xi_tT = xi_tt + P_star_t * (xi_tT - mu_xi - F * xi_tt);
    P_tT  = P_tt + F * (P_tT - P_tp1_t) * F';

    % Store results:
    all_xi_tT(t,:) = xi_tT';
    all_P_tT(t,:)  = reshape(P_tT,1,n_xi^2);
end


