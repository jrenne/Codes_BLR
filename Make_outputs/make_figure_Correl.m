% =========================================================================
% Figure showing consumption growth / inflation correlation
% =========================================================================

indic_k = 5;

T = size(Data_StateSpace.dataset,1);
n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

%     'GDP growth'
n_Y = size(all_xi_tt,2);
vec_dc        = zeros(n_Y,1);
vec_dc(1:n_X) = model_sol.mu_c1;

%     'inflation'
vec_pi  = [model_sol.mu_piX;model_sol.mu_piZ;model_sol.mu_piXX];

Y = all_xi_tt;
[coVar,condCorrel] = compute_condCorrel(model_sol,H,vec_pi,vec_dc,Y);

% Compute gradient of condCorrel
disp("--- Computing conf. interval of conditional correlation ---");
grad_condCorrel = zeros(T,n_Y);
epsilon = .000001;
for i = 1:n_Y
    Y_perturb = Y;
    Y_perturb(:,i) = Y(:,i) + epsilon;
    [coVar,condCorrel_perturb] = ...
        compute_condCorrel(model_sol,H,vec_pi,vec_dc,Y_perturb);
    grad_condCorrel(:,i) = (condCorrel_perturb - condCorrel)/epsilon;
end

kronec = kron(grad_condCorrel,ones(1,n_Y)) .* kron(ones(1,n_Y),grad_condCorrel);
var_condCorrel = (all_P_tt .* kronec) * ones(n_Y^2,1);
std_condCorrel = sqrt(var_condCorrel);

% Create figure
figure('Name','Conditional correlation');

% compute stdv of k (filtering error):
stdv_k = zeros(T,1);
for t = 1:T
    P = reshape(all_P_tt(t,:),n_Y,n_Y);
    stdv_k(t) = sqrt(P(indic_k,indic_k));
end
k = all_xi_tt(:,indic_k);
lower_bound = k - 2*stdv_k;
upper_bound = k + 2*stdv_k;


% Plot factor k
subplot(2,1,1);
plot(dates, k, 'b-', 'LineWidth', 1.5);
title("(a) Factor k(t)",'FontSize', 11);
set(gca, 'FontSize', 9); % increase size of ticks labels
hold on;

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

% Add legend
legend('k factor', '95% confidence band', 'NBER Recession', 'Location', 'best');

% Customize plot
xlabel('Date');
ylabel('Value of factor');
grid on;

% Plot correlation
subplot(2,1,2);
plot(dates, condCorrel, 'b-', 'LineWidth', 1.5);
title("(b) Conditional correlation between consumption growth and inflation",'FontSize', 11);
set(gca, 'FontSize', 9); % increase size of ticks labels
hold on;

lower_bound = condCorrel - 2*std_condCorrel;
upper_bound = condCorrel + 2*std_condCorrel;

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
ylabel('Conditional correlation');
grid on;

% Add legend
legend('Conditional correlation','95% confidence band', 'NBER Recession', 'Location', 'best');

% Save figure in EPS format
figFileName = 'Figures/figure_CondCorrel.eps';
print(figFileName, '-depsc', '-r300');

% Display the plot
hold off;

disp(['Figure saved as ' figFileName]);
