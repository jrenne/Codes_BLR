% =========================================================================
% Figure showing maximum Sharpe ratio
% =========================================================================

Y = all_xi_tt;
maxSR = compute_maxSR(model_sol,Y);

% Compute gradient of maxSR:
disp("--- Computing conf. interval of max Sharpe ratio ---");
grad_condCorrel = zeros(T,n_Y);
epsilon = .000001;
for i = 1:n_Y
    Y_perturb = Y;
    Y_perturb(:,i) = Y(:,i) + epsilon;
    maxSR_perturb = compute_maxSR(model_sol,Y_perturb);
    grad_condCorrel(:,i) = (maxSR_perturb - maxSR)/epsilon;
end

kronec = kron(grad_condCorrel,ones(1,n_Y)) .* kron(ones(1,n_Y),grad_condCorrel);
var_maxSR = (all_P_tt .* kronec) * ones(n_Y^2,1);
std_maxSR = sqrt(var_maxSR);


% Create figure
figure;

% Plot risk aversion
plot(dates, maxSR, 'b-', 'LineWidth', 1.5);
hold on;
ylim([0.5 10])

lower_bound = maxSR - 2*std_maxSR;
upper_bound = maxSR + 2*std_maxSR;

fill([dates' flip(dates)'],...
    [lower_bound' flip(upper_bound)'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% T = size(Data_StateSpace.dataset,1);
% n_X = size(model_sol.PhiQ,1);
% n_Z = size(model_sol.Phi_Z,1);
% n_Y = n_X + n_Z + n_X * (n_X + 1)/2;
% stdv_RA = zeros(T,1);
% for t = 1:T
%     P = reshape(all_P_tt(t,:),n_Y,n_Y);
%     stdv_RA(t) = sqrt(model_sol.mu_gamma1' * P(1:n_X,1:n_X) * model_sol.mu_gamma1);
% end
% lower_bound = riskaversion - 2*stdv_RA;
% upper_bound = riskaversion + 2*stdv_RA;
% %plot(dates,lower_bound);
% %plot(dates,upper_bound);
% fill([dates' flip(dates)'],...
%     [lower_bound' flip(upper_bound)'],...
%     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');


% Shade NBER recession periods
for i = 1:length(recessionStarts)
    recessionPeriod = [recessionStarts(i), recessionEnds(i)];
    xPatch = [recessionPeriod(1), recessionPeriod(2), recessionPeriod(2), recessionPeriod(1)];
    yPatch = [-1000, -1000, 1000, 1000];
    patch(xPatch, yPatch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% Customize plot
xlabel('Date');
ylabel('maximum Sharpe ratio');
grid on;

% Add legend
legend('maximum Sharpe ratio','95% confidence band', 'NBER Recession', 'Location', 'best');


% Save figure in EPS format
figFileName = 'Figures/figure_maxSR.eps';
print(figFileName, '-depsc', '-r300');

% Display the plot
hold off;

disp(['Figure saved as ' figFileName]);
