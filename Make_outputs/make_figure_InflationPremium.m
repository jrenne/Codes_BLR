% =========================================================================
% Figure showing inflation risk premium
% =========================================================================

maturities_in_year = [5;10];
H = max(frequency*maturities_in_year);

% Compute model-implied term premiums:
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


BEIR = nom_yields - real_yields; % Break-Even Inflation Rate
E_Infl = nom_yields_P - real_yields_P; % Expected Inflation
inflation_RP =  BEIR - E_Infl;


figure;

% Plot risk aversion
plot(dates, inflation_RP(:,frequency*maturities_in_year(1)), 'k-', 'LineWidth', 1.5);
hold on;
plot(dates, inflation_RP(:,frequency*maturities_in_year(2)), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);

% Shade NBER recession periods
for i = 1:length(recessionStarts)
    recessionPeriod = [recessionStarts(i), recessionEnds(i)];
    xPatch = [recessionPeriod(1), recessionPeriod(2), recessionPeriod(2), recessionPeriod(1)];
    yPatch = [min(ylim), min(ylim), max(ylim), max(ylim)];
    patch(xPatch, yPatch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

legend('5-year Inflation risk premium','10-year Inflation risk premium',...
    'Location', 'best', 'FontSize', 12);

% Customize plot
xlabel('Date');
ylabel('Inflation risk premium, in percent');
grid on;


% Save the figure in EPS format
epsFileName = 'Figures/figure_InflationPremium.eps';
saveas(gcf, epsFileName, 'epsc');

% Display a message confirming the save
disp(['Figure saved as ' epsFileName]);

