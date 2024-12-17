% =========================================================================
% Figure showing risk aversion coefficient
% =========================================================================

% Prepare the state space model:
[StateSpace,xi_00,P_00,~,~,~,~,~,~,...
    eta0r,eta1r,eta0rn,eta1rn] = prepare_State_Space(model_sol,Data_StateSpace);

% Employ Kalman filter:
n_X = size(model_sol.PhiQ,1);
[all_xi_tt,all_P_tt,logl] = kalman(StateSpace,Data_StateSpace.dataset,xi_00,P_00,n_X);
%[all_xi_tt,all_P_tt] = kalman_smoother(StateSpace,Data_StateSpace.dataset,xi_00,P_00,n_X);

% Plot fit of measured variables:
T = size(Data_StateSpace.dataset,1);
fitted_variables =  ones(T,1) * StateSpace.A + all_xi_tt * StateSpace.H;


figure('Name','Model fit','WindowState','maximized');

count = 0;
for(i = 1:size(Data_StateSpace.dataset,2))
    count = count + 1;
    subplot(4,4,count);
    plot(dates,fitted_variables(:,i), 'k-', 'LineWidth', 2,'Color', [.5 .5 .5]);
    hold on;
    plot(dates,Data_StateSpace.dataset(:,i),'x', 'LineWidth', 2,'Color','black');
    %rmse = sqrt(mean((fitted_variables(:,i)-Data_StateSpace.dataset(:,i)).^2,'omitnan'));
    mae = mean(abs(fitted_variables(:,i)-Data_StateSpace.dataset(:,i)),'omitnan');
    %title([names_of_variables{i} ': ' num2str(mae,'%.2f') ' / ' num2str(sqrt(StateSpace.R(end,i)),'%.2f')],'FontSize', 14);
    title([names_of_variables{i} ' (' num2str(mae,'%.2f') ')'],'FontSize', 14);
    set(gca, 'FontSize', 15); % increase size of ticks labels
    grid on;
    axis tight;
end

% Display the plot
hold off;

if indic_save_output == 1
    % Save figure in EPS format
    figFileName = 'Figures/figure_ModelFit.eps';
    print(figFileName, '-depsc', '-r300');
    disp(['Figure saved as ' figFileName]);
end
