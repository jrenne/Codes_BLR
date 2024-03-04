close all
clear
clc

global FILTER max_abs_param moments frequency Data_StateSpace;


% =========================================================================
% Estimation setup
% =========================================================================
% Re-estimate or used saved results?
indic_estim_moments = 0; % if 1, re-estimate the model
indic_estim_MLE     = 1; % if 1, re-estimate the model
% -------------------------------------------------------------------------
indic_add_moments   = 0;   % minimize -logL + Moment distance
nb_loops_moments    = 1;   % number of estimation loops - moment-fitting approach
nb_loops_MLE        = 20;   % number of estimation loops - MLE approach
nb_iterations_MLE   = 300; % number of iteration for each use of the simplex
indic_save          = 1;   % if == 1, then results are saved
% =========================================================================

% Add paths:
addpath('Data/');
addpath('Plots/');
addpath('Estimation/');
addpath('Procedures/');
addpath('Make_outputs/');

% ========== Load data ====================================================
data_analysis;

% ========== Define initial model =========================================
Make_ini_model;
% Construct and solve model:
model_sol = make_model_sol(model);

% Produce charts:
% plots_simul;


% =========================================================================
% Step 1 --- Moment approach
% =========================================================================

% Define moments to target:
Moments2target;


% ========== Moments fitting ==============================================

max_abs_param = 20; % maximum value of absolute values of parameters

% Define parameters to be optimized on:
FILTER = 0 * model_sol.param + 1;
FILTER([6;8;15;16;17;18]) = 0; % stdv_z, stdv_k, mu_c, mu_pi, mu_gamma0, delta

% Create vector of parameters:
sub_parameters = model.param(FILTER==1);

% Computation of distance at initial parameters:
compute_distance(sub_parameters,model_sol);

f = @(x)compute_distance(x,model_sol);

if indic_estim_moments == 1

    % Set optimization options:
    options = optimset('Display','iter','PlotFcns', ...
        @optimplotfval,'MaxIter',100);

    % Display message box:
    h=msgbox('Please wait. Calculation in progress...');
    for i = 1:nb_loops_moments
        [sub_parameters,fval] = ...
            fminsearch(f,sub_parameters,options);
    end
    % Close message box:
    if isvalid(h); delete(h); end

    % Update model:

    new_parameters = model_sol.param;
    new_parameters(FILTER==1) = sub_parameters;
    model_sol.param = new_parameters;
    model_sol = make_model_sol(model_sol);

    % Print moments:
    indic_disp = 1; % print table
    M = print_moments(sub_parameters,model_sol,indic_disp);

    if indic_save == 1
        save("results/save_moment_approach.mat","model_sol");
    end
else
    load("results/save_moment_approach.mat");
end


% =========================================================================
% Step 2 --- Kalman filter
% =========================================================================

if (indic_estim_moments == 0)&&(indic_estim_MLE == 1)
    % In that case, use Make_ini_model to get initial values
    Make_ini_model;
    % Construct and solve model:
    model_sol = make_model_sol(model);
end

% ========== Kalman filter estimation =====================================
names_of_variables = {'Real consumption of all goods and services',...
    'CPI all ',...
    'output gap',...
    'YIELD3M','YIELD05','YIELD10',...
    'TIPSY05','TIPSY10',...
    'CPI10',...%'RGDP10',...
    'BOND10','BILL10'};
Data_StateSpace = struct;
Data_StateSpace.maturties_nomyields_in_years = [.25;5;10];
Data_StateSpace.maturties_reayields_in_years = [5;10];
Data_StateSpace.names_of_variables = names_of_variables;
% call script:
Make_dataset_observed;


% Maximize log-likelihood: ------------------------------------------------

% Define parameters to be optimized on:
FILTER = 0 * model_sol.param + 1;
FILTER([6;8;9;10;15;16;17;20;21]) = 0;
% stdv_z, stdv_k, mu_gamma1_g, mu_gamma1_z, mu_c, mu_pi, mu_gamma0,
% mu_kappa, sigma_pi
FILTER_MLE = FILTER;

% Create vector of parameters:
sub_parameters =  model_sol.param(FILTER==1);

% Computation of distance at initial parameters:
disp(compute_logl(sub_parameters,Data_StateSpace,model_sol,indic_add_moments));

f = @(x)compute_logl(x,Data_StateSpace,model_sol,indic_add_moments);

if indic_estim_MLE == 1

    % Set optimization options:
    options = optimset('Display','iter','PlotFcns', ...
        @optimplotfval,'MaxIter',nb_iterations_MLE);
    % Display message box:
    %h=msgbox('Please wait.Calculation in progress...');
    for i = 1:nb_loops_MLE
        [sub_parameters,fval] = ...
            fminsearch(f,sub_parameters,options);
    end
    % Close message box:
    %if isvalid(h); delete(h); end

    % Update complete vector of parameters:
    new_parameters = model_sol.param;
    new_parameters(FILTER==1) = sub_parameters;
    model_sol.param = new_parameters;
    model_sol = make_model_sol(model_sol);

    if indic_save == 1
        save("results/save_MLE_approach.mat","model_sol");
    end
else
    load("results/save_MLE_approach.mat");
end

% Computation of distance (to check):
model_sol = make_model_sol(model_sol);
param2print = array2table(round(model_sol.param_transf,4)',...
    'RowNames',model_sol.names_param);
disp(param2print);

FILTER = 0 * model_sol.param + 1; % take all parameters
sub_parameters = model_sol.param(FILTER==1);
disp(compute_logl(sub_parameters,Data_StateSpace,model_sol, ...
    indic_add_moments));

% Print moments:
indic_disp = 1; % print table
M = print_moments(sub_parameters,model_sol,indic_disp);

%plots_simul;

% Plot fit obtained with Kalman filter:
% plot_KF;

% return

%Prepare figures: ---------------------------------------------------------
make_figure_model_fit; % model fit
make_figure_TP; % term premiums
make_figure_Correl; % conditional correlations
make_figure_RA; % Risk aversion
make_figure_RAproxies; % Compare RA with proxies
make_figure_decomp_pi; % decompositions of pi
make_figure_UncondTS; % unconditional term structures of interest
make_figure_InflationPremium; % Inflation risk premiums
make_figure_maxSR; % maximum Sharpe ratio
make_figure_IRF; % Impulse response functions
make_figure_Correl_IRP; % Impulse response functions

%Prepare tables: ----------------------------------------------------------
make_table_moments;
make_table_desc_stat;
make_table_param;
make_table_decompTP;
