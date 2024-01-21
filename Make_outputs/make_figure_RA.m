% =========================================================================
% Figure showing risk aversion coefficient
% =========================================================================

% Compute risk aversion:
riskaversion = model_sol.mu_gamma0*ones(T,1)+all_xi_tt(:,1:n_X)*model_sol.mu_gamma1;

% Create figure
figure;

% Plot risk aversion
plot(dates, riskaversion, 'b-', 'LineWidth', 1.5);
hold on;

T = size(Data_StateSpace.dataset,1);
n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;
stdv_RA = zeros(T,1);
for t = 1:T
    P = reshape(all_P_tt(t,:),n_Y,n_Y);
    stdv_RA(t) = sqrt(model_sol.mu_gamma1' * P(1:n_X,1:n_X) * model_sol.mu_gamma1);
end
lower_bound = riskaversion - 2*stdv_RA;
upper_bound = riskaversion + 2*stdv_RA;
%plot(dates,lower_bound);
%plot(dates,upper_bound);
fill([dates' flip(dates)'],...
    [lower_bound' flip(upper_bound)'],...
    'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');


% Shade NBER recession periods
for i = 1:length(recessionStarts)
    recessionPeriod = [recessionStarts(i), recessionEnds(i)];
    xPatch = [recessionPeriod(1), recessionPeriod(2), recessionPeriod(2), recessionPeriod(1)];
    yPatch = [min(ylim), min(ylim), max(ylim), max(ylim)];
    patch(xPatch, yPatch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% Customize plot
xlabel('Date');
ylabel('Risk aversion');
grid on;

% Add legend
legend('Risk aversion', '95% confidence band', 'NBER Recession', 'Location', 'best');

% Save figure in EPS format
figFileName = 'Figures/figure_RiskAversion.eps';
print(figFileName, '-depsc', '-r300');

% Display the plot
hold off;

disp(['Figure saved as ' figFileName]);
