% =========================================================================
% Figure showing IRFs
% Type 1: Decomposition of nominal rates
% Type 2: Decomposition of real rates
% Type 3: Decomposition of break-even rates
% =========================================================================

matur = 40; % maturity of yield considered
maxH = 10; % time after shock

T = size(Data_StateSpace.dataset,1);
n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;
n_eps = size(model_sol.Sigma,2);

% Pricing -----------------------------------------------------------------

maturities_in_year = 10;
Hmatur = max(frequency*maturities_in_year);

% Compute model-implied term premiums:
[A,B,A4r,B4r] = compute_AB(model_sol,Hmatur);
[An,Bn,Cn,Dn,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,Hmatur);
% Compute specification of rates under Expectation Hypothesis (EH):
[E,V,A1,B1] = compute_EV(model_sol,1); % needed to compute expectations
[APr,BPr,APrn,BPrn] = compute_AB_real_nom_P(model_sol,Hmatur,...
    eta0r,eta1r,eta0rn,eta1rn,A1,B1);

B4r  = [B4r;0*C4rn;0*D4rn];
B4rn = [B4rn;C4rn;D4rn];


% Dynamics ----------------------------------------------------------------

H = 1;
[E,V,A1,B1,Theta0,Theta1,AH,BH] = compute_EV(model_sol,H);
EX = E(1:n_X);
EZ = E((n_X+1):(n_X+n_Z));

vecSigmaZ = model_sol.Gamma0 + model_sol.Gamma1 * EX;
SigmaZ = reshape(vecSigmaZ,n_Z,n_eps);

%% all shocks

% Select shocks:
select_shocks = 1:n_eps; % all

% Number of variables that has changing values:
nb_condition_variable = 6; % for k
% values of k such such that set for kappa = mu_kappa + k = [mu_kappa,0,-mu_kappa]
%values_condit_variable = [0,-model_sol.mu_kappa,-2*model_sol.mu_kappa];
% values of k such such that set for kappa = mu_kappa + k = [-0.5,0,0.5]
values_condit_variable = [-0.75-model_sol.mu_kappa,-model_sol.mu_kappa,0.75-model_sol.mu_kappa];

for type_IRF = 1:3

    switch type_IRF
        case 1
            f1 = figure('Name','IRF real','WindowState','maximized');
            %f1.Position(3:4) = [700 700];
        case 2
            f2 = figure('Name','IRF nominal','WindowState','maximized');
            %f2.Position(3:4) = [700 700];
        case 3
            f3 = figure('Name','IRF BEIR','WindowState','maximized');
            %f3.Position(3:4) = [700 700];
    end

    for type_condit_var = 1:max(size(values_condit_variable))
        
        EXtemp = EX;
        EXtemp(nb_condition_variable) = values_condit_variable(type_condit_var);
        Y0 = [EXtemp;EZ;vech(EXtemp * EXtemp')];

        for ii = 1:max(size(select_shocks))

            n_shock = select_shocks(ii);

            shock = zeros(n_eps,1);
            shock(n_shock) = 1;

            X1 = EXtemp + model_sol.Sigma * shock;
            Z1 = EZ + SigmaZ * shock;
            Y1 = [X1;Z1;vech(X1 * X1')];
            %[Y0 Y1]

            allY0 = zeros(n_Y,maxH+1);
            allY1 = zeros(n_Y,maxH+1);

            all_real_yields_0   = zeros(Hmatur,maxH+1);
            all_nom_yields_0    = zeros(Hmatur,maxH+1);
            all_real_yields_P_0 = zeros(Hmatur,maxH+1);
            all_nom_yields_P_0  = zeros(Hmatur,maxH+1);
            all_real_yields_1   = zeros(Hmatur,maxH+1);
            all_nom_yields_1    = zeros(Hmatur,maxH+1);
            all_real_yields_P_1 = zeros(Hmatur,maxH+1);
            all_nom_yields_P_1  = zeros(Hmatur,maxH+1);

            allY0(:,1) = Y0;
            allY1(:,1) = Y1;

            all_real_yields_0(:,1)   = 100 * frequency * Y0' * B4r;
            all_nom_yields_0(:,1)    = 100 * frequency * Y0' * B4rn;
            all_real_yields_P_0(:,1) = 100 * frequency * Y0' * BPr;
            all_nom_yields_P_0(:,1)  = 100 * frequency * Y0' * BPrn;
            all_real_yields_1(:,1)   = 100 * frequency * Y1' * B4r;
            all_nom_yields_1(:,1)    = 100 * frequency * Y1' * B4rn;
            all_real_yields_P_1(:,1) = 100 * frequency * Y1' * BPr;
            all_nom_yields_P_1(:,1)  = 100 * frequency * Y1' * BPrn;


            Y0h_1 = Y0;
            Y1h_1 = Y1;
            for h = 1:maxH
                Y0h = A1 + B1 * Y0h_1;
                Y1h = A1 + B1 * Y1h_1;

                allY0(:,h+1) = Y0h;
                allY1(:,h+1) = Y1h;

                all_real_yields_0(:,h+1)   = 100 * frequency * Y0h' * B4r;
                all_nom_yields_0(:,h+1)    = 100 * frequency * Y0h' * B4rn;
                all_real_yields_P_0(:,h+1) = 100 * frequency * Y0h' * BPr;
                all_nom_yields_P_0(:,h+1)  = 100 * frequency * Y0h' * BPrn;
                all_real_yields_1(:,h+1)   = 100 * frequency * Y1h' * B4r;
                all_nom_yields_1(:,h+1)    = 100 * frequency * Y1h' * B4rn;
                all_real_yields_P_1(:,h+1) = 100 * frequency * Y1h' * BPr;
                all_nom_yields_P_1(:,h+1)  = 100 * frequency * Y1h' * BPrn;

                Y0h_1 = Y0h;
                Y1h_1 = Y1h;
            end

            IRF = allY1 - allY0;

            %n_var = 1;
            %plot(IRF(n_var,:)');

            IRF_real_yields   = all_real_yields_1 - all_real_yields_0;
            IRF_nom_yields    = all_nom_yields_1  - all_nom_yields_0;
            IRF_real_yields_P = all_real_yields_P_1 - all_real_yields_P_0;
            IRF_nom_yields_P  = all_nom_yields_P_1  - all_nom_yields_P_0;

            switch type_IRF
                case 1
                    yield   = IRF_real_yields(matur,:)';
                    yield_P = IRF_real_yields_P(matur,:)';
                    name_figure = '_real';
                case 2
                    yield   = IRF_nom_yields(matur,:)';
                    yield_P = IRF_nom_yields_P(matur,:)';
                    name_figure = '_nominal';
                case 3
                    yield   = IRF_nom_yields(matur,:)' - IRF_real_yields(matur,:)';
                    yield_P = IRF_nom_yields_P(matur,:)' - IRF_real_yields_P(matur,:)';
                    name_figure = '_BEIR';
            end

            riskpremium = yield - yield_P;

            subplot(max(size(select_shocks)),max(size(values_condit_variable)),...
                (ii-1)*max(size(values_condit_variable))+type_condit_var);
            plot(0:maxH,100*yield, 'k-', 'LineWidth', 1.5);
            hold on;
            plot(0:maxH,100*yield_P,'x', 'LineWidth', 1,'Color', [.5 .5 .5]);
            plot(0:maxH,100*riskpremium,'o', 'LineWidth', 1,'Color','black');
            hold off;

            title(strcat(model_sol.names_shocks_Latex{n_shock},' shock, with $\kappa=', ...
                string(round(values_condit_variable(type_condit_var) + model_sol.mu_kappa,2)),'$'), 'Interpreter', 'latex');

            set(gca, 'FontSize', 12); % increase size of ticks labels

            % Add legend
            %legend('rate','expectations','premiums', 'Location', 'best', 'FontSize', 11);

        end

    end
    
    if indic_save_output == 1
        % Save the figure in EPS format
        epsFileName = ['Figures/figure_IRF_Type' name_figure '.eps'];
        saveas(gcf, epsFileName, 'epsc');
    
        % Display a message confirming the save
        disp(['Figure saved as ' epsFileName]);
    end

end

%% chart for eps_w
select_shocks = 4;

figure('Name','IRF eps_w','WindowState','maximized');

for type_IRF = 1:3

    for type_condit_var = 1:max(size(values_condit_variable))
        
        EXtemp = EX;
        EXtemp(nb_condition_variable) = values_condit_variable(type_condit_var);
        Y0 = [EXtemp;EZ;vech(EXtemp * EXtemp')];

        for ii = 1:max(size(select_shocks))

            n_shock = select_shocks(ii);

            shock = zeros(n_eps,1);
            shock(n_shock) = 1;

            X1 = EXtemp + model_sol.Sigma * shock;
            Z1 = EZ + SigmaZ * shock;
            Y1 = [X1;Z1;vech(X1 * X1')];
            %[Y0 Y1]

            allY0 = zeros(n_Y,maxH+1);
            allY1 = zeros(n_Y,maxH+1);

            all_real_yields_0   = zeros(Hmatur,maxH+1);
            all_nom_yields_0    = zeros(Hmatur,maxH+1);
            all_real_yields_P_0 = zeros(Hmatur,maxH+1);
            all_nom_yields_P_0  = zeros(Hmatur,maxH+1);
            all_real_yields_1   = zeros(Hmatur,maxH+1);
            all_nom_yields_1    = zeros(Hmatur,maxH+1);
            all_real_yields_P_1 = zeros(Hmatur,maxH+1);
            all_nom_yields_P_1  = zeros(Hmatur,maxH+1);

            allY0(:,1) = Y0;
            allY1(:,1) = Y1;

            all_real_yields_0(:,1)   = 100 * frequency * Y0' * B4r;
            all_nom_yields_0(:,1)    = 100 * frequency * Y0' * B4rn;
            all_real_yields_P_0(:,1) = 100 * frequency * Y0' * BPr;
            all_nom_yields_P_0(:,1)  = 100 * frequency * Y0' * BPrn;
            all_real_yields_1(:,1)   = 100 * frequency * Y1' * B4r;
            all_nom_yields_1(:,1)    = 100 * frequency * Y1' * B4rn;
            all_real_yields_P_1(:,1) = 100 * frequency * Y1' * BPr;
            all_nom_yields_P_1(:,1)  = 100 * frequency * Y1' * BPrn;


            Y0h_1 = Y0;
            Y1h_1 = Y1;
            for h = 1:maxH
                Y0h = A1 + B1 * Y0h_1;
                Y1h = A1 + B1 * Y1h_1;

                allY0(:,h+1) = Y0h;
                allY1(:,h+1) = Y1h;

                all_real_yields_0(:,h+1)   = 100 * frequency * Y0h' * B4r;
                all_nom_yields_0(:,h+1)    = 100 * frequency * Y0h' * B4rn;
                all_real_yields_P_0(:,h+1) = 100 * frequency * Y0h' * BPr;
                all_nom_yields_P_0(:,h+1)  = 100 * frequency * Y0h' * BPrn;
                all_real_yields_1(:,h+1)   = 100 * frequency * Y1h' * B4r;
                all_nom_yields_1(:,h+1)    = 100 * frequency * Y1h' * B4rn;
                all_real_yields_P_1(:,h+1) = 100 * frequency * Y1h' * BPr;
                all_nom_yields_P_1(:,h+1)  = 100 * frequency * Y1h' * BPrn;

                Y0h_1 = Y0h;
                Y1h_1 = Y1h;
            end

            IRF = allY1 - allY0;

            %n_var = 1;
            %plot(IRF(n_var,:)');

            IRF_real_yields   = all_real_yields_1 - all_real_yields_0;
            IRF_nom_yields    = all_nom_yields_1  - all_nom_yields_0;
            IRF_real_yields_P = all_real_yields_P_1 - all_real_yields_P_0;
            IRF_nom_yields_P  = all_nom_yields_P_1  - all_nom_yields_P_0;

            switch type_IRF
                case 1
                    yield   = IRF_real_yields(matur,:)';
                    yield_P = IRF_real_yields_P(matur,:)';
                    rate_name = 'Real rate';
                case 2
                    yield   = IRF_nom_yields(matur,:)' - IRF_real_yields(matur,:)';
                    yield_P = IRF_nom_yields_P(matur,:)' - IRF_real_yields_P(matur,:)';
                    rate_name = 'BEIR';
                case 3
                    yield   = IRF_nom_yields(matur,:)';
                    yield_P = IRF_nom_yields_P(matur,:)';
                    rate_name = 'Nominal rate';
            end

            riskpremium = yield - yield_P;

            subplot(3,max(size(values_condit_variable)),(type_IRF-1)*max(size(values_condit_variable))+type_condit_var);
            plot(0:maxH,100*yield, 'k-', 'LineWidth', 2);
            hold on;
            plot(0:maxH,100*yield_P,'kx', 'MarkerSize', 12);
            plot(0:maxH,100*riskpremium,'ko', 'MarkerSize', 12);
            hold off;

            title(strcat(rate_name,', with $\kappa=', ...
                string(round(values_condit_variable(type_condit_var) + model_sol.mu_kappa,2)),'$'), 'Interpreter', 'latex');

            if type_IRF == 2 && type_condit_var == max(size(values_condit_variable))
                % Add legend
                legend('rate','expectations','premiums', 'Location', 'southeast');
            end

            set(gca, 'FontSize', 20); % increase size of ticks labels

        end

    end

end

if indic_save_output == 1
    % Save the figure in EPS format
    epsFileName = 'Figures/figure_IRF_shock_w.eps';
    saveas(gcf, epsFileName, 'epsc');
    
    % Display a message confirming the save
    disp(['Figure saved as ' epsFileName]);
end


%% Slides for eps_w

select_shocks = 4;

for type_IRF = 1:3

    switch type_IRF
        case 1
            figure('Name','IRF real','WindowState','maximized');
        case 2
            figure('Name','IRF nominal','WindowState','maximized');
        case 3
            figure('Name','IRF BEIR','WindowState','maximized');
    end

    for type_condit_var = 1:max(size(values_condit_variable))

        EXtemp = EX;
        EXtemp(nb_condition_variable) = values_condit_variable(type_condit_var);
        Y0 = [EXtemp;EZ;vech(EXtemp * EXtemp')];

        for ii = 1:max(size(select_shocks))

            n_shock = select_shocks(ii);

            shock = zeros(n_eps,1);
            shock(n_shock) = 1;

            X1 = EXtemp + model_sol.Sigma * shock;
            Z1 = EZ + SigmaZ * shock;
            Y1 = [X1;Z1;vech(X1 * X1')];
            %[Y0 Y1]

            allY0 = zeros(n_Y,maxH+1);
            allY1 = zeros(n_Y,maxH+1);

            all_real_yields_0   = zeros(Hmatur,maxH+1);
            all_nom_yields_0    = zeros(Hmatur,maxH+1);
            all_real_yields_P_0 = zeros(Hmatur,maxH+1);
            all_nom_yields_P_0  = zeros(Hmatur,maxH+1);
            all_real_yields_1   = zeros(Hmatur,maxH+1);
            all_nom_yields_1    = zeros(Hmatur,maxH+1);
            all_real_yields_P_1 = zeros(Hmatur,maxH+1);
            all_nom_yields_P_1  = zeros(Hmatur,maxH+1);

            allY0(:,1) = Y0;
            allY1(:,1) = Y1;

            all_real_yields_0(:,1)   = 100 * frequency * Y0' * B4r;
            all_nom_yields_0(:,1)    = 100 * frequency * Y0' * B4rn;
            all_real_yields_P_0(:,1) = 100 * frequency * Y0' * BPr;
            all_nom_yields_P_0(:,1)  = 100 * frequency * Y0' * BPrn;
            all_real_yields_1(:,1)   = 100 * frequency * Y1' * B4r;
            all_nom_yields_1(:,1)    = 100 * frequency * Y1' * B4rn;
            all_real_yields_P_1(:,1) = 100 * frequency * Y1' * BPr;
            all_nom_yields_P_1(:,1)  = 100 * frequency * Y1' * BPrn;


            Y0h_1 = Y0;
            Y1h_1 = Y1;
            for h = 1:maxH
                Y0h = A1 + B1 * Y0h_1;
                Y1h = A1 + B1 * Y1h_1;

                allY0(:,h+1) = Y0h;
                allY1(:,h+1) = Y1h;

                all_real_yields_0(:,h+1)   = 100 * frequency * Y0h' * B4r;
                all_nom_yields_0(:,h+1)    = 100 * frequency * Y0h' * B4rn;
                all_real_yields_P_0(:,h+1) = 100 * frequency * Y0h' * BPr;
                all_nom_yields_P_0(:,h+1)  = 100 * frequency * Y0h' * BPrn;
                all_real_yields_1(:,h+1)   = 100 * frequency * Y1h' * B4r;
                all_nom_yields_1(:,h+1)    = 100 * frequency * Y1h' * B4rn;
                all_real_yields_P_1(:,h+1) = 100 * frequency * Y1h' * BPr;
                all_nom_yields_P_1(:,h+1)  = 100 * frequency * Y1h' * BPrn;

                Y0h_1 = Y0h;
                Y1h_1 = Y1h;
            end

            IRF = allY1 - allY0;

            %n_var = 1;
            %plot(IRF(n_var,:)');

            IRF_real_yields   = all_real_yields_1 - all_real_yields_0;
            IRF_nom_yields    = all_nom_yields_1  - all_nom_yields_0;
            IRF_real_yields_P = all_real_yields_P_1 - all_real_yields_P_0;
            IRF_nom_yields_P  = all_nom_yields_P_1  - all_nom_yields_P_0;

            switch type_IRF
                case 1
                    yield   = IRF_real_yields(matur,:)';
                    yield_P = IRF_real_yields_P(matur,:)';
                    name_figure = '_real';
                case 2
                    yield   = IRF_nom_yields(matur,:)';
                    yield_P = IRF_nom_yields_P(matur,:)';
                    name_figure = '_nominal';
                case 3
                    yield   = IRF_nom_yields(matur,:)' - IRF_real_yields(matur,:)';
                    yield_P = IRF_nom_yields_P(matur,:)' - IRF_real_yields_P(matur,:)';
                    name_figure = '_BEIR';
            end

            riskpremium = yield - yield_P;

            subplot(max(size(select_shocks)),max(size(values_condit_variable)),...
                (ii-1)*max(size(values_condit_variable))+type_condit_var);
            plot(0:maxH,100*yield, 'k-', 'LineWidth', 2);
            hold on;
            plot(0:maxH,100*yield_P,'kx', 'MarkerSize', 15);
            plot(0:maxH,100*riskpremium,'ko', 'MarkerSize', 15);
            hold off;

            title(strcat(model_sol.names_shocks_Latex{n_shock},' shock, with $\kappa=', ...
                string(round(values_condit_variable(type_condit_var) + model_sol.mu_kappa,2)),'$'), 'Interpreter', 'latex');

            if type_condit_var == max(size(values_condit_variable))
                % Add legend
                legend('rate','expectations','premiums', 'Location', 'southeast');
            end

            set(gca, 'FontSize', 30); % increase size of ticks labels

        end

    end

    if indic_save_output == 1
        % Save the figure in EPS format
        epsFileName = ['Figures/figure_IRF_Type' name_figure '_slide.eps'];
        saveas(gcf, epsFileName, 'epsc');
    
        % Display a message confirming the save
        disp(['Figure saved as ' epsFileName]);
    end

end


%% chart for eps_k
select_shocks = 6;

figure('Name','IRF eps_k','WindowState','maximized');

for type_IRF = 1:3

    shock = zeros(n_eps,1);
    shock(select_shocks) = 1;

    X1 = EX + model_sol.Sigma * shock;
    Z1 = EZ + SigmaZ * shock;
    Y1 = [X1;Z1;vech(X1 * X1')];

    all_real_yields_1   = zeros(Hmatur,maxH+1);
    all_nom_yields_1    = zeros(Hmatur,maxH+1);
    all_real_yields_P_1 = zeros(Hmatur,maxH+1);
    all_nom_yields_P_1  = zeros(Hmatur,maxH+1);

    all_real_yields_1(:,1)   = 100 * frequency * Y1' * B4r;
    all_nom_yields_1(:,1)    = 100 * frequency * Y1' * B4rn;
    all_real_yields_P_1(:,1) = 100 * frequency * Y1' * BPr;
    all_nom_yields_P_1(:,1)  = 100 * frequency * Y1' * BPrn;


    Y1h_1 = Y1;
    for h = 1:maxH
        Y1h = A1 + B1 * Y1h_1;

        allY1(:,h+1) = Y1h;

        all_real_yields_1(:,h+1)   = 100 * frequency * Y1h' * B4r;
        all_nom_yields_1(:,h+1)    = 100 * frequency * Y1h' * B4rn;
        all_real_yields_P_1(:,h+1) = 100 * frequency * Y1h' * BPr;
        all_nom_yields_P_1(:,h+1)  = 100 * frequency * Y1h' * BPrn;

        Y1h_1 = Y1h;
    end

    IRF_real_yields   = all_real_yields_1;
    IRF_nom_yields    = all_nom_yields_1;
    IRF_real_yields_P = all_real_yields_P_1;
    IRF_nom_yields_P  = all_nom_yields_P_1;

    switch type_IRF
        case 1
            yield   = IRF_real_yields(matur,:)';
            yield_P = IRF_real_yields_P(matur,:)';
            rate_name = 'Real rate';
        case 2
            yield   = IRF_nom_yields(matur,:)' - IRF_real_yields(matur,:)';
            yield_P = IRF_nom_yields_P(matur,:)' - IRF_real_yields_P(matur,:)';
            rate_name = 'BEIR';
        case 3
            yield   = IRF_nom_yields(matur,:)';
            yield_P = IRF_nom_yields_P(matur,:)';
            rate_name = 'Nominal rate';
    end

    riskpremium = yield - yield_P;

    subplot(2,3,type_IRF);
    plot(0:maxH,100*round(yield,12), 'k-', 'LineWidth', 2);
    hold on;
    plot(0:maxH,100*round(yield_P,12),'kx', 'MarkerSize', 12);
    plot(0:maxH,100*round(riskpremium,12),'ko', 'MarkerSize', 12);
    hold off;
    
    title(rate_name, 'Interpreter', 'latex');

    if type_IRF == 1
        ylim([-1 1]);
        % Add legend
        legend('rate','expectations','premiums', 'Location', 'southwest');
    end

    set(gca, 'FontSize', 20); % increase size of ticks labels

end

if indic_save_output == 1
    % Save the figure in EPS format
    epsFileName = 'Figures/figure_IRF_shock_k.eps';
    saveas(gcf, epsFileName, 'epsc');
    
    % Display a message confirming the save
    disp(['Figure saved as ' epsFileName]);
end
