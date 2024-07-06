% =========================================================================
% Figure comparing RA with proxies
% =========================================================================

% Compute risk aversion:
riskaversion = model_sol.mu_gamma0*ones(T,1) + all_xi_tt(:,1:n_X)*model_sol.mu_gamma1;

% % Load Shiller Excess CAPE: last observation refers to Q4 2019
% X_EXCAPE = xlsread('Data/ShillerCAPE.xlsx','Quarterly');

% Load PSS risk perception: last observation refers to Q4 2023
X_PSS = xlsread('Data/pvs_updated_Jun2023.xlsx','Quarterly');
X_PSS = X_PSS(1:end-16,:); % shorten to Q4 2019

% Load BBM risk apetite: last observation refers to Q4 2023
X_BBM = xlsread('Data/risk_index_bbmdata.xlsx','Quarterly');
X_BBM = X_BBM(1:end-16,:); % shorten to Q4 2019

%AllRA = [riskaversion X_EXCAPE(end-T+1:end,2)];
AllRA = riskaversion;
AllRA(:,end+1:end+2) = nan(T,2);
AllRA(end-length(X_PSS)+1:end,end-1) = X_PSS(:,3);
AllRA(end-length(X_BBM)+1:end,end)   = X_BBM(:,3);

% zscore by Loop due to missings
AllRA_ZSC=NaN(size(AllRA));
for rr=1:size(AllRA,2)
    RAL=AllRA(:,rr);
    AllRA_ZSC(:,rr)=(RAL - nanmean(RAL))/nanstd(RAL);
end

figure;

plot(dates,AllRA_ZSC(:,1), 'k-', 'LineWidth', 2.0);
hold on
plot(dates,AllRA_ZSC(:,2), 'b-.', 'LineWidth', 2.0);  % Grey and thick
plot(dates,AllRA_ZSC(:,3), 'r--', 'LineWidth', 2.0);  % Blue and dashed
%plot(dates,AllRA_ZSC(:,4), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2.0);   % Red and dotted

% Customize plot
xlabel('Date');
ylabel('Series (standardized)');
grid on;
set(gca, 'FontSize', 12); % increase size of ticks labels

% Add legend
%legend('Risk aversion','Shiller Excess CAPE','Location', 'best');
legend('BLR','(negative) risk perception, Pflueger et al.', ...
    '(negative) risk appetite, Bauer et al.','Location', 'northwest', 'FontSize', 11);

hold off;

if indic_save_output == 1
    % Save the figure in EPS format
    epsFileName = 'Figures/figure_RiskAversionComparison.eps';
    saveas(gcf, epsFileName, 'epsc');
    
    % Display a message confirming the save
    disp(['Figure saved as ' epsFileName]);
end

