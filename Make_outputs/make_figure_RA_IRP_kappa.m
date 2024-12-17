% =========================================================================
% Figure showing inflation risk premium vs kappa with risk aversion isolines
% =========================================================================

maturities_in_year = 10;
H = max(frequency*maturities_in_year);

% Compute model-implied term premiums:
[A,B,A4r,B4r] = compute_AB(model_sol,H);
[An,Bn,Cn,Dn,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,H);
% Compute model-implied yields:
real_yields = 100 * frequency * (ones(T,1) * A4r  + X * B4r);
nom_yields  = 100 * frequency * (ones(T,1) * A4rn + X * B4rn + Z * C4rn + XX * D4rn);

% Compute specification of rates under Expectation Hypothesis (EH):
[E,V,A1,B1] = compute_EV(model_sol,1); % needed to compute expectations
[APr,BPr,APrn,BPrn] = compute_AB_real_nom_P(model_sol,H,...
    eta0r,eta1r,eta0rn,eta1rn,A1,B1);
% Compute model-implied yields under EH:
real_yields_P = 100 * frequency * (ones(T,1) * APr  + all_xi_tt * BPr);
nom_yields_P  = 100 * frequency * (ones(T,1) * APrn + all_xi_tt * BPrn);

BEIR   = nom_yields - real_yields; % Break-Even Inflation Rate
E_Infl = nom_yields_P - real_yields_P; % Expected Inflation
inflation_RP =  BEIR - E_Infl;

B4r   = [B4r;0*C4rn;0*D4rn];
B4rn  = [B4rn;C4rn;D4rn];
B4Pr  = BPr;
B4Prn = BPrn;

IRP_constant = 100 * frequency * (A4rn-A4r - (APrn-APr));
IRP_loading  = 100 * frequency * (B4rn-B4r - (BPrn-BPr));

% select states
k_values = -2:0.1:2;
k_values = k_values - model_sol.mu_kappa;
K = length(k_values);
indic_w = 4;
indic_m = 5;
indic_k = 6;
kappa = model_sol.mu_kappa*ones(T,1) + all_xi_tt(:,indic_k);

% Simulate different isolines
Z_simul = zeros(K,n_Z);

% Simulate model for different kappas with gamma_t = mu_gamma for all t
X_simul1 = zeros(K,n_X);
X_simul1(:,indic_k) = k_values';
X_simul1(:,[indic_w,indic_m]) = 0;
XX_simul1 = zeros(K,n_X*(1+n_X)/2);
for k = 1:K
    XX_simul1(k,:) = vech(X_simul1(k,:)' * X_simul1(k,:))';
end
Y_simul1 = [X_simul1 Z_simul XX_simul1];

% Simulate model for different kappas with gamma_t = min(gamma_t) for all t
rho_w   = model_sol.Phi(indic_w,indic_w);
rho_m   = model_sol.Phi(indic_m,indic_m);
sigma_w = model_sol.Sigma(indic_w,indic_w);
sigma_m = model_sol.Sigma(indic_m,indic_m);

var_w = sigma_w^2/(1-rho_w^2);
var_m = sigma_m^2/(1-rho_m^2);

minRA = 1;%min(riskaversion);
X_simul2 = zeros(K,n_X);
X_simul2(:,indic_k) = k_values';
X_simul2(:,indic_w) = (-model_sol.mu_gamma0 + minRA) * var_w/(var_w+var_m);
X_simul2(:,indic_m) = (-model_sol.mu_gamma0 + minRA) * var_m/(var_w+var_m);
XX_simul2 = zeros(K,n_X*(1+n_X)/2);
for k = 1:K
    XX_simul2(k,:) = vech(X_simul2(k,:)' * X_simul2(k,:))';
end
Y_simul2 = [X_simul2 Z_simul XX_simul2];

% Simulate model for different kappas with gamma_t = max(gamma_t) for all t
maxRA = 50;%min(riskaversion);
X_simul3 = zeros(K,n_X);
X_simul3(:,indic_k) = k_values';
X_simul3(:,indic_w) = (-model_sol.mu_gamma0 + maxRA) * var_w/(var_w+var_m);
X_simul3(:,indic_m) = (-model_sol.mu_gamma0 + maxRA) * var_m/(var_w+var_m);
XX_simul3 = zeros(K,n_X*(1+n_X)/2);
for k = 1:K
    XX_simul3(k,:) = vech(X_simul3(k,:)' * X_simul3(k,:))';
end
Y_simul3 = [X_simul3 Z_simul XX_simul3];

% Calculate kappa, IRP and risk aversion for all simulations
kappa_simul1 = model_sol.mu_kappa*ones(K,1) + Y_simul1(:,indic_k);
kappa_simul2 = model_sol.mu_kappa*ones(K,1) + Y_simul2(:,indic_k);
kappa_simul3 = model_sol.mu_kappa*ones(K,1) + Y_simul3(:,indic_k);

IRP_simul1 = ones(K,1)*IRP_constant + Y_simul1*IRP_loading;
IRP_simul2 = ones(K,1)*IRP_constant + Y_simul2*IRP_loading;
IRP_simul3 = ones(K,1)*IRP_constant + Y_simul3*IRP_loading;

riskaversion1 = model_sol.mu_gamma0*ones(K,1) + Y_simul1(:,1:n_X)*model_sol.mu_gamma1;
riskaversion2 = model_sol.mu_gamma0*ones(K,1) + Y_simul2(:,1:n_X)*model_sol.mu_gamma1;
riskaversion3 = model_sol.mu_gamma0*ones(K,1) + Y_simul3(:,1:n_X)*model_sol.mu_gamma1;

% Plot scatterplot with isolines
figure;
scatter(kappa,inflation_RP(:,frequency*maturities_in_year), 50,'b');
hold on;
plot(kappa_simul2,IRP_simul2(:,frequency*maturities_in_year),'k:', 'LineWidth', 2);
plot(kappa_simul3,IRP_simul3(:,frequency*maturities_in_year),'k--', 'LineWidth', 2);
plot(kappa_simul1,IRP_simul1(:,frequency*maturities_in_year),'k-', 'LineWidth', 2);
hold off;
%legend('model estimates', 'isoline at \gamma = min(\gamma_t)','isoline at \gamma = max(\gamma_t)','isoline at \gamma = \mu_\gamma', 'Location', 'northeast', 'FontSize', 11);
legend('model estimates', 'isoline at \gamma = 1','isoline at \gamma = 50','isoline at \gamma = \mu_\gamma', 'Location', 'northeast', 'FontSize', 11);
ylabel('10-year inflation risk premium, in percent');
xlabel('kappa');
grid on;
xlim([-1 0.5]);
ylim([-3 5]);
set(gca, 'FontSize', 12); % increase size of ticks labels
box on;

if indic_save_output == 1
    % Save the figure in EPS format
    epsFileName = 'Figures/figure_RA_IRP_kappa.eps';
    saveas(gcf, epsFileName, 'epsc');

    % Display a message confirming the save
    disp(['Figure saved as ' epsFileName]);
end
