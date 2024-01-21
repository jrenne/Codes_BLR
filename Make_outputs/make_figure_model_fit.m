% =========================================================================
% Figure showing risk aversion coefficient
% =========================================================================

% Prepare the state space model:
[StateSpace,xi_00,P_00,~,~,~,~,~,~,...
    eta0r,eta1r,eta0rn,eta1rn] = prepare_State_Space(model_sol,Data_StateSpace);

% Employ Kalman filter:
[all_xi_tt,all_P_tt,logl] = kalman(StateSpace,Data_StateSpace.dataset,xi_00,P_00);
%[all_xi_tt,all_P_tt] = kalman_smoother(StateSpace,Data_StateSpace.dataset,xi_00,P_00);

% Plot fit of measured variables:
T = size(Data_StateSpace.dataset,1);
fitted_variables =  ones(T,1) * StateSpace.A + all_xi_tt * StateSpace.H;


figure('Name','Model fit','WindowState','maximized');

count = 0;
for(i = 1:size(Data_StateSpace.dataset,2))
    count = count + 1;
    subplot(6,2,count);
    plot(dates,fitted_variables(:,i), 'k-', 'LineWidth', 2,'Color', [.5 .5 .5]);
    hold on;
    plot(dates,Data_StateSpace.dataset(:,i),'x', 'LineWidth', 2,'Color','black');
    title(names_of_variables{i},'FontSize', 16);
    set(gca, 'FontSize', 13); % increase size of ticks labels
    grid on;
end

% Save figure in EPS format
figFileName = 'Figures/figure_ModelFit.eps';
print(figFileName, '-depsc', '-r300');

% Display the plot
hold off;

disp(['Figure saved as ' figFileName]);
