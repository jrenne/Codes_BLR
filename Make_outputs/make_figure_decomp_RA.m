% =========================================================================
% Figure showing risk aversion decomposition
% =========================================================================

maturities_in_year = [1;10];

H = max(frequency*maturities_in_year);

% Compute model-implied term premiums:
[A,B,A4r,B4r] = compute_AB(model_sol,H);

% Compute model-implied yields:
real_yields = 100 * frequency * (ones(T,1)*A4r + X * B4r);

% Compute specification of rates under Expectation Hypothesis (EH):
[E,V,A1,B1] = compute_EV(model_sol,1); % needed to compute expectations
[APr,BPr,APrn,BPrn] = compute_AB_real_nom_P(model_sol,H,...
    eta0r,eta1r,eta0rn,eta1rn,A1,B1);
% Compute model-implied yields under EH:
real_yields_P = 100 * frequency * (ones(T,1)*APr + all_xi_tt * BPr);


% Real term premium
real_TP = real_yields - real_yields_P;

% Compute risk aversion:
riskaversion = model_sol.mu_gamma0*ones(T,1) + all_xi_tt(:,1:n_X)*model_sol.mu_gamma1;

% Create figure
figure('WindowState','maximized');

m = 1;
subplot(1,3,1);
plot(dates, riskaversion, 'b-', 'LineWidth', 1.5);
yyaxis right
hold on;
plot(dates, real_TP(:,maturities_in_year(m)*frequency));
hold off;
legend('risk aversion (lhs)', 'real term premium (rhs)')
title(['Risk aversion vs ' num2str(maturities_in_year(m)) '-year real TP']);

m = 2;
subplot(1,3,2);
plot(dates, riskaversion, 'b-', 'LineWidth', 1.5);
yyaxis right
hold on;
plot(dates, real_TP(:,maturities_in_year(m)*frequency));
hold off;
legend('risk aversion (lhs)', 'real term premium (rhs)')
title(['Risk aversion vs ' num2str(maturities_in_year(m)) '-year real TP']);

subplot(1,3,3);
bar(dates,[model_sol.mu_gamma0*ones(T,1) all_xi_tt(:,[4 5])],'stacked');
hold on;
plot(dates, riskaversion, 'b-', 'LineWidth', 3.0);
hold off;
legend('mean','high-frequency','low-frequency','risk aversion');
title('Decomposition of risk aversion');
