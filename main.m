close all
clear
clc

global FILTER max_abs_param moments frequency Data_StateSpace;


%% ========================================================================
% Estimation setup
% =========================================================================
% Re-estimate or used saved results?
indic_estim_moments = 0; % if 1, re-estimate the model
indic_estim_MLE     = 0; % if 1, re-estimate the model
% -------------------------------------------------------------------------
indic_add_moments   = 0;   % minimize -logL + Moment distance
nb_loops_moments    = 3;   % number of estimation loops - moment-fitting approach
nb_loops_MLE        = 30;  % number of estimation loops - MLE approach
nb_iterations_MLE   = 100; % number of iteration for each use of the simplex
indic_save_model    = 0;   % if == 1, then mode results are saved
indic_save_output   = 0;   % if == 1, then figures & tables are saved
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

%% ========================================================================
% Step 1 --- Moment approach
% =========================================================================

% Define moments to target:
Moments2target;


% ========== Moments fitting ==============================================

max_abs_param = 20; % maximum value of absolute values of parameters

% Define parameters to be optimized on:
FILTER = 0 * model_sol.param + 1;
FILTER([9,10,16,17]) = 0; %sigma_m,sigma_k,mu_c,mu_gamma

% Create vector of parameters:
sub_parameters = model.param(FILTER==1);

% Computation of distance at initial parameters:
%compute_distance(sub_parameters,model_sol);

f = @(x)compute_distance(x,model_sol);

if indic_estim_moments == 1

    % Set optimization options:
    options = optimset('Display','iter','PlotFcns', ...
        @optimplotfval,'MaxIter',100);

    for i = 1:nb_loops_moments
        [sub_parameters,fval] = ...
            fminsearch(f,sub_parameters,options);
    end

    % Update model:
    new_parameters = model_sol.param;
    new_parameters(FILTER==1) = sub_parameters;
    model_sol.param = new_parameters;
    model_sol = make_model_sol(model_sol);

    % Print moments:
    indic_disp = 1; % print table
    M = print_moments(sub_parameters,model_sol,indic_disp);

    if indic_save_model == 1
        save("results/save_moment_approach.mat","model_sol");
    end
else
    load("results/save_moment_approach.mat");
end

%% ========================================================================
% Step 2 --- Kalman filter
% =========================================================================

if (indic_estim_moments == 0)&&(indic_estim_MLE ~= 0)
    % In that case, use Make_ini_model to get initial values
    Make_ini_model;
    % Construct and solve model:
    model_sol = make_model_sol(model);
end

% ========== Kalman filter estimation =====================================
names_of_variables = {'Consumption growth',...
    'CPI inflation', 'Output gap',...
    'YLD3M','YLD02','YLD05','YLD10','YLD20', ...
    'TIPSY02','TIPSY05','TIPSY10','TIPSY20', ...
    'RTR','PTR','BILL10','CPI10'};
%    'YIELD3M','YIELD02','YIELD05','YIELD10','YIELD20', ... 
%    'TIPSY02','TIPSY05','REALR10','REALR20', ... % 'TIPSY10','TIPSY20',
%    'TIPSY02','REAL05','REAL10','REAL20', ...
%    'PTR','CPI10','RTR','SurveyRTP10'};
%    'RGDP10','BOND10','BILL10'};
Data_StateSpace = struct;
Data_StateSpace.maturties_nomyields_in_years = [.25;2;5;10;20];
Data_StateSpace.maturties_reayields_in_years = [2;5;10;20];
Data_StateSpace.names_of_variables = names_of_variables;
% call script:
Make_dataset_observed;

% Maximize log-likelihood: ------------------------------------------------

% Define parameters to be optimized on:
FILTER = 0 * model_sol.param + 1;
FILTER([9,10,16,17]) = 0; %sigma_m,sigma_k,mu_c,mu_gamma
FILTER_MLE = FILTER;

% Create vector of parameters:
sub_parameters =  model_sol.param(FILTER==1);

% Computation of distance at initial parameters:
%disp(compute_logl(sub_parameters,Data_StateSpace,model_sol,indic_add_moments));

f = @(x)compute_logl(x,Data_StateSpace,model_sol,indic_add_moments);

if indic_estim_MLE == 1

    % Set optimization options:
    options = optimset('Display','iter','PlotFcns', @optimplotfval,'MaxIter',nb_iterations_MLE);

    % First, run simplex optimiser (runs faster)
    for i = 1:nb_loops_MLE
        disp(i);
        [sub_parameters,fval] = fminsearch(f,sub_parameters,options);
    end

    % Second, iterate through simplex and gradient-based optimiser
    for i = 1:nb_loops_MLE
        disp(i);
        [sub_parameters,~] = fminsearch(f,sub_parameters,options);
        [sub_parameters,fval,~,~,~,Hessian] = fminunc(f,sub_parameters,options);
    end
    model_sol.Hessian = Hessian;

    % min_flag    = 0;
    % count_loops = 0;
    % while min_flag ~= 1
    %     count_loops = count_loops + 1;
    %     [sub_parameters,~] = fminsearch(f,sub_parameters,options);
    %     [sub_parameters,fval,~,~,~,Hessian] = fminunc(f,sub_parameters,options);
    %     model_sol.Hessian = Hessian;
    %     try
    %         chol(Hessian);
    %         min_flag = 1;
    %         disp(['The Hessian is positive definite. Required extra loops: ' num2str(count_loops)]);
    %         return;
    %     catch
    %         disp(['The Hessian was not positive definite. Tried extra loops: ' num2str(count_loops)]);
    %     end
    % end

    % Update complete vector of parameters:
    new_parameters = model_sol.param;
    new_parameters(FILTER==1) = sub_parameters;
    model_sol.param = new_parameters;
    model_sol = make_model_sol(model_sol);

    if indic_save_model == 1
        save("results/save_MLE_approach.mat","model_sol");
    end
else
    load("results/save_MLE_approach.mat");
end

%% Results

% Computation of distance (to check):
% model_sol = make_model_sol(model_sol);
% param2print = array2table(round(model_sol.param_transf,4)',...
%     'RowNames',model_sol.names_param);
% disp(param2print);

FILTER = 0 * model_sol.param + 1; % take all parameters
sub_parameters = model_sol.param(FILTER==1);
disp(f(sub_parameters));

% Print moments:
indic_disp = 1; % print table
M = print_moments(sub_parameters,model_sol,indic_disp);

%plots_simul;

% Plot fit obtained with Kalman filter:
%plot_KF;

% Calculate covariance-matrix of parameters
compute_stdev_params;

%Prepare figures: ---------------------------------------------------------
make_figure_model_fit;     % model fit
make_figure_TP;            % nominal term premiums
make_figure_Correl;        % conditional correlations
make_figure_RA;            % Risk aversion
make_figure_RAproxies;     % compare RA with proxies
make_figure_decomp_RA;     % decomposition of risk aversion
make_figure_decomp_pi;     % decompositions of pi
make_figure_UncondTS;      % unconditional term structures of interest
make_figure_Premiums;      % Inflation risk and real term premiums
make_figure_decomp_rates;  % decomposition of yields, tips and beirs
make_figure_IRF;           % Impulse response functions
make_figure_Correl_IRP;    % Correlation vs IRP
make_figure_RA_IRP_kappa;  % Scatterplot between kappa and IRP
make_figure_RA_RTP;        % Risk aversion vs RTP
% make_figure_slope;         % hysteresis effect on real slope
% make_figure_maxSR;         % maximum Sharpe ratio
% make_figure_Correl_IRP_3D; % Correlation vs IRP

%Prepare tables: ----------------------------------------------------------
make_table_moments;
make_table_desc_stat;
make_table_param;
make_table_decompTP;
