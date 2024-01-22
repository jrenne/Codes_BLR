% =========================================================================
% Figure showing Term Premiums
% =========================================================================

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

X  = all_xi_tt(:,1:n_X);
Z  = all_xi_tt(:,(n_X+1):(n_X+n_Z));
XX = all_xi_tt(:,(n_X+n_Z+1):end);


[X_TP,TXT_TP,~] = xlsread('Data/TermPremiaUSForCheck.xlsx','TP_Comparison_US');
TP_sel = [4,5,8,9]; % TP_05y_ACM, TP_05y_KW, TP_10y_ACM, TP_10y_KW
TP_label = TXT_TP(1,TP_sel-1);
TP = nan(T,length(TP_sel));
delay = 8; %the TP data set starts 8 quarters later
TP(delay+1:end,:) = X_TP(1:T-delay,TP_sel); 

indic_variable = find(strcmp(data_names,"YIELD10"));
YIELD10 = data(:,indic_variable);
indic_variable = find(strcmp(data_names,"YIELD05"));
YIELD05 = data(:,indic_variable);

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

figure('Name','TP comparison', ...
    'WindowState','maximized');

subplot(2,2,1);
plot(dates,nom_yields(:,frequency*maturities_in_year(1)), 'k-', 'LineWidth', 2);
hold on;
% Plot the second line with grey color and dotted style
plot(dates, nom_yields_P(:, frequency * maturities_in_year(1)), ...
    '--', 'LineWidth', 2,'Color', [.5 .5 .5]);
plot(dates,YIELD05,'x', 'LineWidth', 2,'Color','black');
title('5-year yields',...%'Interpreter', 'latex', 'FontSize', 20);
    'FontSize', 20);
% Add legend with increased FontSize
legend('Fitted yield', 'Model-implied expectation component','Data',...
    'Location', 'northeast', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 16); % increase size of ticks labels

subplot(2,2,2);
plot(dates,nom_yields(:,frequency*maturities_in_year(1))-...
    nom_yields_P(:,frequency*maturities_in_year(1)), 'b-', 'LineWidth', 2);  % Blue and thick
hold on;
plot(dates,TP(:,1), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);  % Grey and thick
plot(dates,TP(:,2), 'r--', 'LineWidth', 2);  % Red and thick
legend('BLR','KW','ACM','Location', 'northeast', 'FontSize', 16);
title('5-year term premiums', 'FontSize', 20);
grid on
set(gca, 'FontSize', 16); % increase size of ticks labels

subplot(2,2,3);
plot(dates,nom_yields(:,frequency*maturities_in_year(2)), 'k-', 'LineWidth', 2);
hold on;
% Plot the second line with grey color and dotted style
plot(dates, nom_yields_P(:, frequency * maturities_in_year(2)), ...
    '--', 'LineWidth', 2,'Color', [.5 .5 .5]);
plot(dates,YIELD10,'x', 'LineWidth', 2,'Color','black');
title('10-year yields', 'FontSize', 20);
% Add legend with increased FontSize
legend('Fitted yield', 'Model-implied expectation component','Data',...
    'Location', 'northeast', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 16); % increase size of ticks labels


subplot(2,2,4);
plot(dates,nom_yields(:,frequency*maturities_in_year(2))-...
    nom_yields_P(:,frequency*maturities_in_year(2)), 'b-', 'LineWidth', 2);  % Blue and thick
hold on;
plot(dates,TP(:,3), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);  % Grey and thick
plot(dates,TP(:,4), 'r--', 'LineWidth', 2);  % Red and thick
legend('BLR','KW','ACM','Location', 'northeast', 'FontSize', 16);
title('10-year term premiums', 'FontSize', 20);
grid on;
set(gca, 'FontSize', 16); % increase size of ticks labels

hold off;

% Save the figure in EPS format
epsFileName = 'Figures/figure_TP.eps';
saveas(gcf, epsFileName, 'epsc');

% Display a message confirming the save
disp(['Figure saved as ' epsFileName]);

