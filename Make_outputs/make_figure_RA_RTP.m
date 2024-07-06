% =========================================================================
% Figure showing real term premia & risk aversion
% =========================================================================

% Compute risk aversion:
riskaversion = model_sol.mu_gamma0*ones(T,1) + all_xi_tt(:,1:n_X)*model_sol.mu_gamma1;


maturities_in_year = [5;10];
H = max(frequency*maturities_in_year);

% Compute model-implied term premiums:
[A,B,A4r,B4r] = compute_AB(model_sol,H);
[An,Bn,Cn,Dn,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,H);

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

figure;

yyaxis left
plot(dates, 100 * real_TP(:,frequency*maturities_in_year(2)), 'b-', 'LineWidth', 2);
hold on;
plot(dates, 100 * real_TP(:,frequency*maturities_in_year(1)), 'b:','LineWidth', 2);
ylabel('Real term premium, in bps','Color','b');
yticks(-50:50:200);
ylim([min(yticks),max(yticks)]);

yyaxis right
plot(dates, riskaversion,'r', 'LineWidth', 2);
ylabel('Risk aversion','Color','r');
yticks(-20:20:80);
ylim([min(yticks),max(yticks)]);
hold off;

set(gca, 'FontSize', 12); % increase size of ticks labels

% % Shade NBER recession periods
% for i = 1:length(recessionStarts)
%     recessionPeriod = [recessionStarts(i), recessionEnds(i)];
%     xPatch = [recessionPeriod(1), recessionPeriod(2), recessionPeriod(2), recessionPeriod(1)];
%     yPatch = [min(ylim), min(ylim), max(ylim), max(ylim)];
%     patch(xPatch, yPatch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% end

legend('10-year RTP','5-year RTP','Risk aversion','Location', 'southeast', 'FontSize', 11);

% Customize plot
xlabel('Date');
grid on;

if indic_save_output == 1
    % Save the figure in EPS format
    epsFileName = 'Figures/figure_RA_RTP.eps';
    saveas(gcf, epsFileName, 'epsc');
    
    % Display a message confirming the save
    disp(['Figure saved as ' epsFileName]);
end
