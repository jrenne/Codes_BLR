% =========================================================================
% Figure showing consumption growth / inflation correlation
% =========================================================================

horiz_condCorrel = 20;

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

%     'output gap'
vec_z    = zeros(n_Y,1);
vec_z(2) = 1;
%vec_z(3) = -1;

%     'cyclical inflation'
vec_pi_tilde  = zeros(n_Y,1);
vec_pi_tilde(n_X+2) = 1;

Y = all_xi_tt;

model_sol_new = model_sol;
%model_sol_new.param(2) = -20;
%model_sol_new.param(4) = -20;
%model_sol_new.param(8) = -20;
model_sol_new = make_model_sol(model_sol_new);

condCorrel       = nan(T,horiz_condCorrel);
condCorrel_cycle = nan(T,horiz_condCorrel);
coVar            = nan(T,horiz_condCorrel);
coVar_cycle      = nan(T,horiz_condCorrel);
for h = 1:horiz_condCorrel
    [coVar_temp,condCorrel_temp] = compute_condCorrel(model_sol_new,h,vec_pi,vec_dc,Y);
    condCorrel(:,h) = condCorrel_temp;
    coVar(:,h) = coVar_temp;
    [coVar_temp,condCorrel_temp] = compute_condCorrel(model_sol_new,h,vec_pi_tilde,vec_z,Y);
    condCorrel_cycle(:,h) = condCorrel_temp;
    coVar_cycle(:,h) = coVar_temp;
end

% Create figure
figure('Name','Conditional correlation 3D');

subplot(1,2,1);
% Plot correlation
%plot(dates, condCorrel, 'b-', 'LineWidth', 1.5);
mesh(1:horiz_condCorrel,dates,condCorrel);
title("Conditional correlation between consumption growth and inflation",'FontSize', 11);
set(gca, 'FontSize', 9); % increase size of ticks labels
set(gca, 'yDir','reverse')
% Customize plot
xlabel('Forward horizon (quarters)');
ylabel('Date');
zlabel('Conditional correlation');
grid on;

subplot(1,2,2);
% Plot correlation
%plot(dates, condCorrel, 'b-', 'LineWidth', 1.5);
mesh(1:horiz_condCorrel,dates,condCorrel_cycle);
title("Conditional correlation between cyclical components (z & pi tilde)",'FontSize', 11);
set(gca, 'FontSize', 9); % increase size of ticks labels
set(gca, 'yDir','reverse')
% Customize plot
xlabel('Forward horizon (quarters)');
ylabel('Date');
zlabel('Conditional correlation');
grid on;

% Save figure in EPS format
figFileName = 'Figures/figure_CondCorrel_3D.eps';
print(figFileName, '-depsc', '-r300');

disp(['Figure saved as ' figFileName]);
