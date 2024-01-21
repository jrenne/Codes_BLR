% =========================================================================
% Figure showing risk unconditional term structures of interest rates
% =========================================================================


% model_sol.names_param;
% model_sol.param(19) = 0.9578;
% model_sol = make_model_sol(model_sol);


% Compute moments of Y:
H = 1;
[E,V,A1,B1,~,~,~,BH,~,~] = compute_EV(model_sol,H);

% Compute uncond yield curves:
max_H = 10;
H = max_H * frequency;
[A,B,A4r,B4r] = compute_AB(model_sol,H);
[An,Bn,Cn,Dn,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,H);

loadings4r  = frequency * [B4r;0*C4rn;0*D4rn];
loadings4rn = frequency * [B4rn;C4rn;D4rn];

% compute stdv of yields:
all_stdv_r  = zeros(H,1);
all_stdv_rn = zeros(H,1);
for h = 1:H
    all_stdv_r(h)  = sqrt(loadings4r(:,h)' * V * loadings4r(:,h));
    all_stdv_rn(h) = sqrt(loadings4rn(:,h)' * V * loadings4rn(:,h));
end

uncond_r  = frequency * (A4r'+loadings4r'*E);
uncond_rn = frequency * (A4rn'+loadings4rn'*E);

% Plot average TS of real and nominal yields:
n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

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
