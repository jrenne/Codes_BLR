

% Prepare the state space model:
[StateSpace,xi_00,P_00,~,~,~,~,~,~,...
    eta0r,eta1r,eta0rn,eta1rn] = prepare_State_Space(model_sol,Data_StateSpace);

% Employ Kalman filter:
[all_xi_tt,all_P_tt,logl] = kalman(StateSpace,Data_StateSpace.dataset,xi_00,P_00);
%[all_xi_tt,all_P_tt] = kalman_smoother(StateSpace,Data_StateSpace.dataset,xi_00,P_00);

% Plot fit of measured variables:
T = size(Data_StateSpace.dataset,1);
n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

fitted_variables =  ones(T,1) * StateSpace.A + all_xi_tt * StateSpace.H;

name_figure_fit   = 'Fit of observed variables (post MLE)';
figure('Name',name_figure_fit,'WindowState','maximized');

count = 0;
for(i = 1:size(Data_StateSpace.dataset,2))
    count = count + 1;
    subplot(4,4,count);
    plot(dates,[fitted_variables(:,i) Data_StateSpace.dataset(:,i)]);
    title(names_of_variables{i});
    if names_of_variables{i} == "CPI all "
        title("inflation (trend in blue, cycle in green)");
        hold on
        plot(dates,100 * all_xi_tt(:,n_X+1),'-b');
        plot(dates,100 * all_xi_tt(:,n_X+2),'--g');
        hold off
    end
    if (names_of_variables{i} == "BOND10")|...
            (names_of_variables{i} == "BILL10")|...
            (names_of_variables{i} == "RGDP10")
        hold on
        plot(dates,Data_StateSpace.dataset(:,i),"o");
        hold off
    end
end



% =========================================================================
% Model-implied objects
% =========================================================================

name_figure_addit = 'Model outputs (post MLE)';
figure('Name',name_figure_addit,'WindowState','maximized');


% Risk aversion:

count = 1;
subplot(3,2,count);
riskaversion = model_sol.mu_gamma0*ones(T,1)+all_xi_tt(:,1:n_X)*model_sol.mu_gamma1;
plot(dates,[riskaversion all_xi_tt(:,4)])
title('Risk aversion (w in red)')

% Term premiums

% Compute model-implied real and nominal yield curves:
H = 10 * frequency;

X  = all_xi_tt(:,1:n_X);
Z  = all_xi_tt(:,(n_X+1):(n_X+n_Z));
XX = all_xi_tt(:,(n_X+n_Z+1):end);

[A,B,A4r,B4r] = compute_AB(model_sol,H);
[An,Bn,Cn,Dn,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,H);
% Compute model-implied yields:
real_yields = 100 * frequency * (ones(T,1)*A4r + X * B4r);
nom_yields  = 100 * frequency * (ones(T,1)*A4rn + X * B4rn + Z * C4rn + XX * D4rn);

% Compute specification of rates under Expectation Hypothesis (EH):
[E,V,A1,B1] = compute_EV(model_sol,1); % needed to compute expectations
[APr,BPr,APrn,BPrn] = compute_AB_real_nom_P(model_sol,H,...
    eta0r,eta1r,eta0rn,eta1rn,A1,B1);
% Compute model-implied yields under EH:
real_yields_P = 100 * frequency * (ones(T,1)*APr + all_xi_tt * BPr);
nom_yields_P  = 100 * frequency * (ones(T,1)*APrn + all_xi_tt * BPrn);

maturities_in_year = [5;10];
count = count + 1;
subplot(3,2,count)
plot(dates,[nom_yields(:,frequency*maturities_in_year)-nom_yields_P(:,frequency*maturities_in_year)])
title("10 year nominal term premium (5yr in blue, 10yr in red)")
count = count + 1;
subplot(3,2,count)
plot(dates,[real_yields(:,frequency*maturities_in_year)-real_yields_P(:,frequency*maturities_in_year)])
title("10 year real term premium (5yr in blue, 10yr in red)")


% Inflation

count = count + 1;
subplot(3,2,count)
k = all_xi_tt(:,5);
fk = 1/2 * (1 + k);
gk = 1/2 * (1 - k);
plot(dates,[fk.^2 gk.^2 k]);
title('f(k)^2 in blue (variance share, demand), g(k)^2 in blue [k in yellow]')

% Computation of conditionnal correlation between consumption growth and
%  inflation

H = 4;
[E,V,A,B,Theta0,Theta1,AH,BH,Theta0H,Theta1H] = compute_EV(model_sol,H);

% Compute conditional covariance matrices:
condCov = ones(T,1) * Theta0H' + all_xi_tt * Theta1H';

%     'GDP growth'
n_Y = size(all_xi_tt,2);
vec_dc        = zeros(n_Y,1);
vec_dc(1:n_X) = model_sol.mu_c1;
%     'output gap'
vec_z    = zeros(n_Y,1);
vec_z(2) = 1;
%     'inflation'
vec_pi  = [model_sol.mu_piX;model_sol.mu_piZ;model_sol.mu_piXX];

coVar = condCov .* (ones(T,1) * kron(vec_dc',vec_pi'));
coVar = coVar * ones(n_Y^2,1);

var_dc = condCov .* (ones(T,1) * kron(vec_dc',vec_dc'));
var_dc = var_dc * ones(n_Y^2,1);

var_pi = condCov .* (ones(T,1) * kron(vec_pi',vec_pi'));
var_pi = var_pi * ones(n_Y^2,1);

condCorrel = coVar ./ (sqrt(var_pi) .* sqrt(var_dc));

count = count + 1;
subplot(3,2,count)
plot(condCorrel);
title("4-quarter-ahead correlation between inflation and consump. growth");




