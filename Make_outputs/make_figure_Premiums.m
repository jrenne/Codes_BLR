% =========================================================================
% Figures showing inflation risk and real term premiums
% =========================================================================

% Load DKW estimates (5y and 10y IRP)
X_DKW = xlsread('Data/DKW_updates.xlsx','quarterly');

% Last observation refers to Q4 2019
DKW_IRP = X_DKW(end-T+1:end,2:3);
DKW_LP  = X_DKW(end-T+1:end,4:5);
DKW_RTP = X_DKW(end-T+1:end,6:7);

% % Load AACMY estimates (10y and 5y5y IRP)
% X_AACMY = xlsread('Data/IRP_AACMY.xlsx','quarterly');
% 
% % Last observation refers to Q4 2019
% AACMY_IRP = nan(T,2);
% AACMY_IRP(T-length(X_AACMY)+1:end,:) = X_AACMY(:,2:3);

% Calculate own estimates

maturities_in_year = 10;
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

%BEIR = nom_yields - real_yields; % Break-Even Inflation Rate

%IRP
E_Infl = nom_yields_P - real_yields_P; % Expected Inflation
inflation_RP =  nom_yields - real_yields - E_Infl;

% Real term premium
real_TP = real_yields - real_yields_P;

% Survey-implied premia
indx_SurveyIRP10 = find(strcmp(data_names,{'SurveyIRP10'}));
SurveyIRP10 = data(:,indx_SurveyIRP10);
indx_SurveyRTP10 = find(strcmp(data_names,{'SurveyRTP10'}));
SurveyRTP10 = data(:,indx_SurveyRTP10);

%% Inflation risk premium

figure;
plot(dates, inflation_RP(:,frequency*maturities_in_year(1)), 'k-', 'LineWidth', 2);
hold on;
%plot(dates, inflation_RP(:,frequency*maturities_in_year(2)), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
plot(dates, DKW_IRP(:,2), 'b-.', 'LineWidth', 2);
%plot(dates, AACMY_IRP(:,1), 'LineWidth', 2);
plot(dates, SurveyIRP10, 'r--', 'LineWidth', 2);

set(gca, 'FontSize', 12); % increase size of ticks labels

% Shade NBER recession periods
for i = 1:length(recessionStarts)
    recessionPeriod = [recessionStarts(i), recessionEnds(i)];
    xPatch = [recessionPeriod(1), recessionPeriod(2), recessionPeriod(2), recessionPeriod(1)];
    yPatch = [min(ylim), min(ylim), max(ylim), max(ylim)];
    patch(xPatch, yPatch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% legend('BLR 5y','BLR 10y',...
%     'DKW 5y', 'DKW 10y',...
%     'AACMY 10y',...
%     'Survey-implied 10y obs', ...
%     'Location', 'best', 'FontSize', 11);

legend('BLR','DKW','Survey-implied','Location', 'best', 'FontSize', 11)

% Customize plot
xlabel('Date');
ylabel('Inflation risk premium, in percent');
grid on;

if indic_save_output == 1
    % Save the figure in EPS format
    epsFileName = 'Figures/figure_InflationPremiums.eps';
    saveas(gcf, epsFileName, 'epsc');
    
    % Display a message confirming the save
    disp(['Figure saved as ' epsFileName]);
end

%% Real term premiums

figure;
plot(dates, real_TP(:,frequency*maturities_in_year(1)), 'k-', 'LineWidth', 2);
hold on;
%plot(dates, real_TP(:,frequency*maturities_in_year(2)), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
plot(dates, DKW_RTP(:,2), 'b-.', 'LineWidth', 2);
plot(dates, SurveyRTP10,  'r--','LineWidth', 2);

set(gca, 'FontSize', 12); % increase size of ticks labels

% Shade NBER recession periods
for i = 1:length(recessionStarts)
    recessionPeriod = [recessionStarts(i), recessionEnds(i)];
    xPatch = [recessionPeriod(1), recessionPeriod(2), recessionPeriod(2), recessionPeriod(1)];
    yPatch = [min(ylim), min(ylim), max(ylim), max(ylim)];
    patch(xPatch, yPatch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% legend('BLR 5y','BLR 10y',...
%     'DKW 5y', 'DKW 10y',...
%     'Survey-implied 10y obs', ...
%     'Location', 'best', 'FontSize', 11);

legend('BLR','DKW','Survey-implied','Location', 'best', 'FontSize', 11);

% Customize plot
xlabel('Date');
ylabel('Real term premium, in percent');
grid on;

if indic_save_output == 1
    % Save the figure in EPS format
    epsFileName = 'Figures/figure_RealPremiums.eps';
    saveas(gcf, epsFileName, 'epsc');
    
    % Display a message confirming the save
    disp(['Figure saved as ' epsFileName]);
end
