% =========================================================================
% Figure showing risk unconditional term structures of interest rates
% =========================================================================

% Compute uncondtional yield curves:
H = frequency * 20;
[uncond_r,uncond_rn,uncond_TPr,uncond_TPrn,uncond_IRP,...
    all_stdv_r,all_stdv_rn,all_stdv_TPr,all_stdv_TPrn,all_stdv_IRP] =...
    compute_uncond_yds_TP(model_sol,H);

% Create figure
figure;

plot((1:H)/frequency,uncond_r, 'b-', 'LineWidth', 1.5);
hold on;
plot((1:H)/frequency,uncond_rn, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 12); % increase size of ticks labels

% % Shaded areas with standard deviations
% upper_bound = uncond_r + all_stdv_r;
% lower_bound = uncond_r - all_stdv_r;
% xPatch = [(1:H)/frequency flip((1:H)/frequency)];
% yPatch = [lower_bound' flip(upper_bound)'];
% patch(xPatch, yPatch, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% 
% upper_bound = uncond_rn + all_stdv_rn;
% lower_bound = uncond_rn - all_stdv_rn;
% xPatch = [(1:H)/frequency flip((1:H)/frequency)];
% yPatch = [lower_bound' flip(upper_bound)'];
% patch(xPatch, yPatch, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

ylabel('Interest rates, in percent');
xlabel('Maturity, in years');

slope_data_r = mean(Data_StateSpace.dataset(:,9:12),'omitnan')';
slope_data_n = mean(Data_StateSpace.dataset(:,4:8),'omitnan')';

slope_fit_r = fitted_variables(:,9:12);
slope_fit_r(isnan(Data_StateSpace.dataset(:,9:12))) = nan;
slope_fit_r = mean(slope_fit_r,'omitnan')';

slope_fit_n = fitted_variables(:,4:8);
slope_fit_n(isnan(Data_StateSpace.dataset(:,4:8))) = nan;
slope_fit_n = mean(slope_fit_n,'omitnan')';

grid_n = Data_StateSpace.maturties_nomyields_in_years;
grid_r = Data_StateSpace.maturties_reayields_in_years;

scatter(grid_r,slope_data_r,'b','filled');
scatter(grid_r,slope_fit_r,100,'b');
scatter(grid_n,slope_data_n,'r','filled');
scatter(grid_n,slope_fit_n,100,'r');

legend('Real rates', 'Nominal rates',...%'+/- 1 std dev','+/- 1 std dev',...
    'sample average (real)', 'fitted average (real)', ...
    'sample average (nominal)', 'fitted average (nominal)', ...
    'Location', 'southeast','FontSize', 11);

ylim([-5 10]);

% Display the plot
hold off;

if indic_save_output == 1
    % Save figure in EPS format
    figFileName = 'Figures/figure_UncondTS.eps';
    print(figFileName, '-depsc', '-r300');
    disp(['Figure saved as ' figFileName]);
end