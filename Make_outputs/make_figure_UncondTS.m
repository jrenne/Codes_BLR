% =========================================================================
% Figure showing risk unconditional term structures of interest rates
% =========================================================================

% Compute uncondtional yield curves:
H = frequency * 10;
[uncond_r,uncond_rn,uncond_TPr,uncond_TPrn,uncond_IRP,...
    all_stdv_r,all_stdv_rn,all_stdv_TPr,all_stdv_TPrn,all_stdv_IRP] =...
    compute_uncond_yds_TP(model_sol,H);

% Check what happens when rho_gz = 0 (no hyteresis)
indic_variable = find(strcmp(model_sol.names_param,'rho_gz'));
model_noHysteresis = model_sol;
model_noHysteresis.param(indic_variable) = 0.01;
model_noHysteresis_sol = make_model_sol(model_noHysteresis);
[uncond_r_noHyteresis] = compute_uncond_yds_TP(model_noHysteresis_sol,H);


% Create figure
figure;

plot((1:H)/frequency,uncond_r, 'b-', 'LineWidth', 1.5);
hold on;
plot((1:H)/frequency,uncond_rn, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 12); % increase size of ticks labels

% Shaded areas with standard deviations
upper_bound = uncond_r + all_stdv_r;
lower_bound = uncond_r - all_stdv_r;
xPatch = [(1:H)/frequency flip((1:H)/frequency)];
yPatch = [lower_bound' flip(upper_bound)'];
patch(xPatch, yPatch, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

upper_bound = uncond_rn + all_stdv_rn;
lower_bound = uncond_rn - all_stdv_rn;
xPatch = [(1:H)/frequency flip((1:H)/frequency)];
yPatch = [lower_bound' flip(upper_bound)'];
patch(xPatch, yPatch, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

ylabel('Interest rates, in percent', 'FontSize', 14);
xlabel('Maturity, in years', 'FontSize', 14);

legend('Real rates', 'Nominal rates',...
    '+/- 1 std dev','+/- 1 std dev',...
    'Location', 'northwest', 'FontSize', 12);


% Save figure in EPS format
figFileName = 'Figures/figure_UncondTS.eps';
print(figFileName, '-depsc', '-r300');

% Display the plot
hold off;

disp(['Figure saved as ' figFileName]);
