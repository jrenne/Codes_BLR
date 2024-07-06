close all
clear
clc

global FILTER max_abs_param moments frequency Data_StateSpace;

% Add paths:
addpath('Estimation/');
addpath('Procedures/');
addpath('Make_outputs/');


H = 10;
vec_rho_gz = [.0000 .003];

range = [-100 50];

nb_sigma_z  = 30;
nb_sigma_g  = 30;
min_sigma_g = .00001;
max_sigma_g = .005;
min_sigma_z = .00001;
max_sigma_z = .05;

sigma_g_dot = .00;
sigma_z_dot = .04;

z_values = linspace(min_sigma_z,max_sigma_z,nb_sigma_z);
g_values = linspace(min_sigma_g,max_sigma_g,nb_sigma_g);

% Construct initial model:
Make_ini_model;

model.param_transf.rho_g  = .9;
model.param_transf.rho_z  = .8;

model.param_transf.mu_c    = .02;

model.param_transf.mu_gamma0   = 10;
model.param_transf.mu_gamma1_g = 0;
model.param_transf.mu_gamma1_z = 0;
model.param_transf.sigma_w     = 0;

RA_ini = model.param_transf.mu_gamma0;

levels = linspace(range(1),range(2),40);

f = figure('Name','Real YCslope');
f.Position(3:4) = [1000 700];

customColormap = b2r(-6,3);

for iii = 1:2

    model.param_transf.rho_gz = vec_rho_gz(iii);

    % Initialize an array to store results
    results = zeros(nb_sigma_z,nb_sigma_g);

    % Use parfor loop for parallel computing
    for i = 1:nb_sigma_z

        model.param_transf.sigma_z = z_values(i);

        for j = 1:nb_sigma_g

            model.param_transf.sigma_g = g_values(j);

            model.param   = make_param(model);
            model_sol_new = make_model_sol(model);
            [A,B,A4r,B4r] = compute_AB(model_sol_new,H);

            results(i,j) = A4r(H) - A4r(1);
        end
    end

    % Display the final results
    %disp('Final Results:');
    %disp(results);

    subplot(2,2,iii);

    % Create a custom blue-red colormap
    contourf(g_values,z_values,10000*results,levels);
    caxis(range);

    % Apply the custom colormap
    colormap(customColormap);

    c = colorbar;
    set(gca, 'FontSize', 14);

    % Specify the size of the point
    hold on;
    pointSize = 100;
    scatter(sigma_g_dot, sigma_z_dot, ...
        pointSize, 'ro', 'filled');
    hold off;

    xlabel('$\sigma_g$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$\sigma_z$', 'Interpreter', 'latex', 'FontSize', 16);
    title(['10-year real term premium, with $\rho_{gz}=$ '...
        num2str(model.param_transf.rho_gz)], 'Interpreter', 'latex', 'FontSize', 16);

    if(iii == 1)
        colorbar('off');
    end

    % Prepare IRF
    Phi = model_sol_new.Phi;
    shock = zeros(size(Phi,1),1);
    shock(2) = 1;
    x = shock;

    maxH = 40;
    all_x = zeros(size(Phi,1),maxH+1);
    all_dc = zeros(1,maxH+1);
    all_x(:,1) = x;
    all_dc(1) = model_sol_new.mu_c1' * x;
    for h = 1:maxH
        x = Phi * x;
        all_x(:,h+1) = x;
        all_dc(h+1) = model_sol_new.mu_c1' * x;
    end

    subplot(2,2,2+iii)
    plot(0:maxH,cumsum(all_dc), 'k-', 'LineWidth', 2);
    set(gca, 'FontSize', 14);

    xlabel('time after shock', 'FontSize', 14);
    ylabel('Consumption effect', 'FontSize', 14);
    title(['Consumption effect of output gap shock $\rho_{gz}=$ '...
        num2str(model.param_transf.rho_gz)], 'Interpreter', 'latex', 'FontSize', 16);

end


% Save the figure in EPS format
epsFileName = 'Figures/figure_slope.eps';
saveas(gcf, epsFileName, 'epsc');

% Display a message confirming the save
disp(['Figure saved as ' epsFileName]);




% =========================================================================
% Figure showing inflence of RA on term premiums

model.param_transf.sigma_z = sigma_z_dot;
model.param_transf.sigma_g = sigma_g_dot;

nb_RA_values = 30;
RA_values = linspace(0.5,20,nb_RA_values);

results = zeros(nb_RA_values,2);

f = figure('Name','TP as function of RA');
%f.Position(3:4) = [400 300];
%f.Position(4) = [200];

for iii = 1:2

    model.param_transf.rho_gz = vec_rho_gz(iii);

    % Use parfor loop for parallel computing
    for j = 1:nb_RA_values

        model.param_transf.mu_gamma0 = RA_values(j);

        model.param   = make_param(model);
        model_sol_new = make_model_sol(model);
        [A,B,A4r,B4r] = compute_AB(model_sol_new,H);

        results(j,iii) = A4r(H) - A4r(1);

        % Create model where P = Q:
        model_P = model_sol_new;
        model_P.muQ = 0 * model_P.muQ;
        model_P.PhiQ = model_P.Phi;
        [A_P,B_P,A4r_P,B4r_P] = compute_AB(model_P,H);
        %results(j,iii) = A4r(H) - A4r_P(H);

    end
end

plot(RA_values,10000*results(:,1), 'k--', 'LineWidth', 2);
hold on;
plot(RA_values,10000*results(:,2), 'k-', 'LineWidth', 2);

set(gca, 'FontSize', 14);

ylabel('10-year real term premium (in bps)', 'FontSize', 13);
xlabel('$\gamma$', 'Interpreter', 'latex', 'FontSize', 18);

plot([RA_ini, RA_ini], ylim, 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
hold off;

% Add legend
legend(['without hysteresis effect $(\rho_{gz} =' num2str(vec_rho_gz(1)) ')$'], ...
    ['with hysteresis effect $(\rho_{gz} =' num2str(vec_rho_gz(2)) ')$'],...
    'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 15);

if indic_save_output == 1
    % Save the figure in EPS format
    epsFileName = 'Figures/figure_slope_RA.eps';
    saveas(gcf, epsFileName, 'epsc');
    
    % Display a message confirming the save
    disp(['Figure saved as ' epsFileName]);
end


