function[all_xi_tt,all_P_tt,logl,logl_t,all_Gains] = kalman(StateSpace,DataObs,xi_00,P_00,option)

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

mu_xi = StateSpace.mu_xi;
F     = StateSpace.F;
%R     = StateSpace.R;
A     = StateSpace.A;
H     = StateSpace.H;

n_xi = size(F,1);
n_y  = size(H,2);

T = size(DataObs,1); % number of periods

% Initialize matrices:
all_xi_tt = zeros(T,n_xi);
all_P_tt  = zeros(T,n_xi^2);

all_Gains = NaN(T,n_xi*n_y);

% Initialize vectors (date 0):
xi_t_1t_1 = xi_00;
P_t_1t_1  = P_00;

% Log-likelihood:
%logl   = 0;
logl_t = zeros(T,1);

for t = 1:T
    
    R = diag(StateSpace.R(t,:));

    y_t = DataObs(t,:)';

    Q = functionQ(StateSpace,xi_t_1t_1);

    % Forecasting step:
    xi_tt_1 = mu_xi + F * xi_t_1t_1 ;
    y_tt_1  = A' + H' * xi_tt_1;
    P_tt_1  = F * P_t_1t_1 * F' + Q;

    % Updating step:
    lambda_t = y_t - y_tt_1;
    indic_notNaN = ~isnan(lambda_t);
    if sum(indic_notNaN)==0
        % There is no observation on this date.
        xi_tt    = xi_tt_1;
        P_tt     = P_tt_1;
        Gain_not_reduced = NaN(n_xi,n_y);
    else
        % There is at least one observation on this date.
        Hstar = H(:,indic_notNaN);
        Rstar = R(:,indic_notNaN);
        Rstar = Rstar(indic_notNaN,:);
        nstar = sum(indic_notNaN);
        lambda_tstar = lambda_t(indic_notNaN);
        Gain     = P_tt_1 * Hstar * (Hstar' * P_tt_1 * Hstar + Rstar)^(-1);
        xi_tt    = xi_tt_1 + Gain * lambda_tstar;
        P_tt     = P_tt_1 - Gain * Hstar' * P_tt_1;

        Gain_not_reduced = NaN(n_xi,n_y);
        Gain_not_reduced(:,indic_notNaN) = Gain;
    end

    % Impose constraints on xi components:
    if(~isnan(option))
        xi_tt = impose_constraint(xi_tt,option);
    end

    % Store results:
    all_xi_tt(t,:) = xi_tt';
    all_P_tt(t,:)  = reshape(P_tt,1,n_xi^2);
    all_Gains(t,:) = reshape(Gain_not_reduced,1,n_xi*n_y);

    % Update of log-likelihood:
    % logl = logl - nstar/2 * log(2*pi) - 1/2 * det(Hstar' * P_tt_1 * Hstar + Rstar) - ...
    %     1/2 * lambda_tstar' * (Hstar' * P_tt_1 * Hstar + Rstar)^(-1) * lambda_tstar;

    logl_t(t) = - nstar/2 * log(2*pi) - 1/2 * det(Hstar' * P_tt_1 * Hstar + Rstar) - ...
        1/2 * lambda_tstar' * (Hstar' * P_tt_1 * Hstar + Rstar)^(-1) * lambda_tstar;

    xi_t_1t_1 = xi_tt;
    P_t_1t_1  = P_tt;
end

logl = sum(logl_t);

