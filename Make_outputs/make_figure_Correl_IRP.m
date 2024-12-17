% =========================================================================
% Figure showing inflation risk premium & consumption-inflation correlation
% =========================================================================

maturities_in_year = [5;10];
H = max(frequency*maturities_in_year);
T = size(Data_StateSpace.dataset,1);

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

% Compute Inflation / Growth correlation:
n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);

% Compute conditional correlation:
horiz_condCorrel = 1;
Y = all_xi_tt;

%     'GDP growth'
n_Y = size(all_xi_tt,2);
vec_dc        = zeros(n_Y,1);
vec_dc(1:n_X) = model_sol.mu_c1;

%     'inflation'
vec_pi  = [model_sol.mu_piX;model_sol.mu_piZ;model_sol.mu_piXX];

[~,condCorrel]     = compute_condCorrel(model_sol,horiz_condCorrel,vec_dc,vec_pi,Y);


figure;

yyaxis left
plot(dates, 100 * inflation_RP(:,frequency*maturities_in_year(2)), 'b-', 'LineWidth', 2);
hold on;
plot(dates, 100 * inflation_RP(:,frequency*maturities_in_year(1)), 'b:','LineWidth', 2);
ylabel('Inflation risk premium, in bps','Color','b');
yticks(-400:200:400);
ylim([min(yticks),max(yticks)]);

yyaxis right
plot(dates, condCorrel,'r', 'LineWidth', 2);
ylabel('Conditional correlation','Color','r');
yticks(-1:0.5:1);
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

legend('10-year IRP','5-year IRP','Conditional correlation','Location', 'northeast', 'FontSize', 11);

% Customize plot
xlabel('Date');
grid on;

if indic_save_output == 1
    % Save the figure in EPS format
    epsFileName = 'Figures/figure_Correl_IRP.eps';
    saveas(gcf, epsFileName, 'epsc');
    
    % Display a message confirming the save
    disp(['Figure saved as ' epsFileName]);
end
