% =========================================================================
% Figure showing pi decomposition
% =========================================================================

n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

figure;

indic_variable = find(strcmp(data_names,'CPI all '));
inflation = data(:,indic_variable);
%    'output gap',...

plot(dates,inflation, 'k-', 'LineWidth', 1.5);
hold on
plot(dates,100 * all_xi_tt(:,n_X+1),'-', 'LineWidth', 2,'Color', [.7 .7 .7]);
plot(dates,100 * all_xi_tt(:,n_X+2),'--b', 'LineWidth', 1.5);
hold off;

grid on;

% Add legend
legend('Inflation', 'Exogenous component', 'Cyclical (z-related) component', 'Location', 'best');


% Save the figure in EPS format
epsFileName = 'Figures/figure_decomp_pi.eps';
saveas(gcf, epsFileName, 'epsc');

% Display a message confirming the save
disp(['Figure saved as ' epsFileName]);

