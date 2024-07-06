% =========================================================================
% Figure showing pi decomposition
% =========================================================================

n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

figure;

% indic_variable = strcmp(data_names,'CPI inflation');
% inflation = data(:,indic_variable);

indic_variable = strcmp(Data_StateSpace.names_of_variables,'CPI inflation');
inflation_fit = fitted_variables(:,indic_variable);

plot(dates,inflation_fit, 'k-', 'LineWidth', 1.5);
hold on
plot(dates,100 * (model_sol.mu_pi0+all_xi_tt(:,n_X+1)),'-', 'LineWidth', 2,'Color', [.7 .7 .7]);
plot(dates,100 * all_xi_tt(:,n_X+2),'--b', 'LineWidth', 1.5);
hold off;

grid on;
set(gca, 'FontSize', 12); % increase size of ticks labels

% Add legend
legend('Fitted inflation', 'Model-implied mean + pi star', 'Cyclical (z-related) component', 'Location', 'best', 'FontSize', 11);

if indic_save_output == 1
    % Save the figure in EPS format
    epsFileName = 'Figures/figure_decomp_pi.eps';
    saveas(gcf, epsFileName, 'epsc');
    
    % Display a message confirming the save
    disp(['Figure saved as ' epsFileName]);
end
