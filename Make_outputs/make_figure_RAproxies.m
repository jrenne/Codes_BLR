% =========================================================================
% Figure comparing RA with proxies
% =========================================================================

% WL adds co-plot with other RA proxies
RAProxies = textread('Data/RAProxies.txt');

AllRA=[riskaversion RAProxies];
% zscore by Loop due to missings
AllRA_ZSC=NaN(size(AllRA));
for rr=1:3;
    RAL=AllRA(:,rr);
    AllRA_ZSC(:,rr)=(RAL - nanmean(RAL))/nanstd(RAL);
end;

figure;
%plot(dates,AllRA_ZSC(:,1:3));

plot(dates,AllRA_ZSC(:,1), 'k-', 'LineWidth', 2);
hold on
plot(dates,AllRA_ZSC(:,2), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);  % Grey and thick
plot(dates,AllRA_ZSC(:,3), 'b--', 'LineWidth', 1.5);  % Red and thick

% Customize plot
xlabel('Date');
%ylabel('Risk aversion / Excess CAPE (standardized)');
ylabel('Series (standardized)');
grid off;

% Add legend
%legend('Risk aversion','Shiller Excess CAPE','Location', 'best');
legend('BLR','Shiller Excess CAPE','(negative) risk appetite, Pflueger et al','Location', 'northwest');

hold off;


% Save the figure in EPS format
epsFileName = 'Figures/figure_RiskAversionComparison.eps';
saveas(gcf, epsFileName, 'epsc');

% Display a message confirming the save
disp(['Figure saved as ' epsFileName]);


