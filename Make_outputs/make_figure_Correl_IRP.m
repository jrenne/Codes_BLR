% =========================================================================
% Figure showing inflation risk premium
% =========================================================================

maturities_in_year = [5;10];
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

%B4IRP = (B4rn - B4r) - (B4Prn - B4Pr);
%IRP = 100 * frequency * all_xi_tt * diag(B4IRP(:,40));
%plot(IRP(:,1));
%plot(IRP * ones(22,1));

% Compute Inflation / Growth correlation:

indic_k = 5;

T = size(Data_StateSpace.dataset,1);
n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;
% Compute conditional correlation:
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
vec_z([2,3]) = 1;

%     'inflation'
vec_pi  = [model_sol.mu_piX;model_sol.mu_piZ;model_sol.mu_piXX];

coVar = condCov .* (ones(T,1) * kron(vec_dc',vec_pi'));
coVar = coVar * ones(n_Y^2,1);

var_dc = condCov .* (ones(T,1) * kron(vec_dc',vec_dc'));
var_dc = var_dc * ones(n_Y^2,1);

var_pi = condCov .* (ones(T,1) * kron(vec_pi',vec_pi'));
var_pi = var_pi * ones(n_Y^2,1);

condCorrel = 100*coVar ./ (sqrt(var_pi) .* sqrt(var_dc));


figure;

yyaxis left
plot(dates, 100 * inflation_RP(:,frequency*maturities_in_year(1)), 'b-', 'LineWidth', 1.5);
ylabel('Inflation risk premium, in bps','Color','b');

yyaxis right
plot(dates, condCorrel,'r-', 'LineWidth', 1.5);
ylabel('Conditional correlation, in percent','Color','r');
ylim([-12,12]);

% Plot risk aversion
%hold on;
%hold off;


% Shade NBER recession periods
for i = 1:length(recessionStarts)
    recessionPeriod = [recessionStarts(i), recessionEnds(i)];
    xPatch = [recessionPeriod(1), recessionPeriod(2), recessionPeriod(2), recessionPeriod(1)];
    yPatch = [min(ylim), min(ylim), max(ylim), max(ylim)];
    patch(xPatch, yPatch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

legend('5-year Inflation risk premium','Conditional correlation',...
    'Location', 'best', 'FontSize', 12);

% Customize plot
xlabel('Date');
grid on;


% Save the figure in EPS format
epsFileName = 'Figures/figure_Correl_IRP.eps';
saveas(gcf, epsFileName, 'epsc');

% Display a message confirming the save
disp(['Figure saved as ' epsFileName]);

