% =========================================================================
% Table showing decomposition of term premiums variance
% =========================================================================

Format = '%0.2f';

maturities_in_year = [2;5;10];
H = max(frequency*maturities_in_year);

param2set2zero = {'sigma_w','sigma_k','mu_kappa','rho_gz','sigma_g','sigma_z'};

% Baseline case:
[uncond_r,uncond_rn,uncond_TPr,uncond_TPrn,uncond_IRP,...
    all_stdv_r,all_stdv_rn,all_stdv_TPr,all_stdv_TPrn,all_stdv_IRP] =...
    compute_uncond_yds_TP(model_sol,H);

% Create LaTeX table

nb_maturity = max(size(maturities_in_year));
latexTable = sprintf('\\begin{tabular}{cc');
for i = 1:(3*(1+nb_maturity))
    latexTable = [latexTable 'c'];
end
latexTable = [latexTable '} \hline '];

latexTable = [latexTable '&&\multicolumn{' num2str(nb_maturity) '}{c}{Nominal term prem.}'...
    '&&\multicolumn{' num2str(nb_maturity) '}{c}{Real term prem.}'...
    '&&\multicolumn{' num2str(nb_maturity) '}{c}{Inflation risk prem.}'...
    ];
latexTable = [latexTable '\\ \hline '];

str_maturities = [num2str(maturities_in_year(1)) ' yrs'];
for i = 2:nb_maturity
    str_maturities = [str_maturities '& ' num2str(maturities_in_year(i)) ' yrs'];
end

latexTable = [latexTable '&&' str_maturities '&&' str_maturities '&&' str_maturities '\\ \hline '];

latexTable = [latexTable '\\ && \multicolumn{' num2str(3*nb_maturity+2) '}{c}{\bf{A. Unconditional means of term premiums}}\\ \hline '];

% Baseline term premiums:
latexTable = [latexTable 'baseline &'];
% Nominal:
for i = 1:nb_maturity
    latexTable = [latexTable ' & ' sprintf(Format,...
        uncond_TPrn(frequency * maturities_in_year(i)))];
end
latexTable = [latexTable ' & '];
% Real:
for i = 1:nb_maturity
    latexTable = [latexTable ' & ' sprintf(Format,...
        uncond_TPr(frequency * maturities_in_year(i)))];
end
latexTable = [latexTable ' & '];
% Inflation:
for i = 1:nb_maturity
    latexTable = [latexTable ' & ' sprintf(Format,...
        uncond_IRP(frequency * maturities_in_year(i)))];
end

latexTable = [latexTable ' \\ '];

for(iii = 1:max(size(param2set2zero)))

    % first entry of line:
    indic_param = find(strcmp(model_sol.names_param,param2set2zero(iii)));
    name_latex = model_sol.names_param_Latex{indic_param};
    %latexTable = [latexTable '$' erase(name_latex,'$') '=0$ &' ];

    % Modify model and compute TPs:
    model_sol_new = model_sol;
  
    if strcmp(param2set2zero(iii),'sigma_w')
        model_sol_new.param(indic_param) = -Inf;
        indic_param = find(strcmp(model_sol.names_param,'sigma_m'));
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'constant RA &' ];
    end

    if strcmp(param2set2zero(iii),'sigma_k')
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'constant $corr(\pi,c)$ &' ];
    end

    if strcmp(param2set2zero(iii),'mu_kappa')
        model_sol_new.param(indic_param) = 0;
        indic_param = find(strcmp(model_sol.names_param,'sigma_k'));
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'zero $corr(\pi,c)$ &' ];
    end

    if strcmp(param2set2zero(iii),'rho_gz')
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'no hysteresis &' ];
    end

    if strcmp(param2set2zero(iii),'sigma_g')
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'no trend growth &' ];
    end

    if strcmp(param2set2zero(iii),'sigma_z')
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'no cyclical growth &' ];
    end

    % avoid recalculation of endogenous parameters
    endo_flag = 0;
    model_sol_new = make_model_sol(model_sol_new,endo_flag);

    % Compute term premiums:
    [uncond_r_new,uncond_rn_new,uncond_TPr_new,uncond_TPrn_new,uncond_IRP_new,...
        all_stdv_r_new,all_stdv_rn_new,all_stdv_TPr_new,all_stdv_TPrn_new,all_stdv_IRP_new] =...
        compute_uncond_yds_TP(model_sol_new,H);

    % Nominal:
    for i = 1:nb_maturity
        latexTable = [latexTable ' & ' sprintf(Format,...
            uncond_TPrn_new(frequency * maturities_in_year(i)))];
    end
    latexTable = [latexTable ' & '];
    % Real:
    for i = 1:nb_maturity
        latexTable = [latexTable ' & ' sprintf(Format,...
            uncond_TPr_new(frequency * maturities_in_year(i)))];
    end
    latexTable = [latexTable ' & '];
    % Inflation:
    for i = 1:nb_maturity
        latexTable = [latexTable ' & ' sprintf(Format,...
            uncond_IRP_new(frequency * maturities_in_year(i)))];
    end

    latexTable = [latexTable ' \\ '];

end

latexTable = [latexTable ' \hline '];

latexTable = [latexTable '\\ && \multicolumn{' num2str(3*nb_maturity+2) '}{c}{\bf{B. Unconditional standard deviation of term premiums}}\\ \hline '];

% Baseline term premiums:
latexTable = [latexTable 'baseline &'];
% Nominal:
for i = 1:nb_maturity
    latexTable = [latexTable ' & ' sprintf(Format,...
        all_stdv_TPrn(frequency * maturities_in_year(i)))];
end
latexTable = [latexTable ' & '];
% Real:
for i = 1:nb_maturity
    latexTable = [latexTable ' & ' sprintf(Format,...
        all_stdv_TPr(frequency * maturities_in_year(i)))];
end
latexTable = [latexTable ' & '];
% Inflation:
for i = 1:nb_maturity
    latexTable = [latexTable ' & ' sprintf(Format,...
        all_stdv_IRP(frequency * maturities_in_year(i)))];
end

latexTable = [latexTable ' \\ '];


for(iii = 1:max(size(param2set2zero)))

    % first entry of line:
    indic_param = find(strcmp(model_sol.names_param,param2set2zero(iii)));
    name_latex = model_sol.names_param_Latex{indic_param};
    %latexTable = [latexTable '$' erase(name_latex,'$') '=0$ &' ];

    % Modify model and compute TPs:
    model_sol_new = model_sol;
    
    if strcmp(param2set2zero(iii),'sigma_w')
        model_sol_new.param(indic_param) = -Inf;
        indic_param = find(strcmp(model_sol.names_param,'sigma_m'));
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'constant RA &' ];
    end

    if strcmp(param2set2zero(iii),'sigma_k')
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'constant $corr(\pi,c)$ &' ];
    end

    if strcmp(param2set2zero(iii),'mu_kappa')
        model_sol_new.param(indic_param) = 0;
        indic_param = find(strcmp(model_sol.names_param,'sigma_k'));
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'zero $corr(\pi,c)$ &' ];
    end

    if strcmp(param2set2zero(iii),'rho_gz')
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'no hysteresis&' ];
    end

    if strcmp(param2set2zero(iii),'sigma_g')
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'no trend growth &' ];
    end


    if strcmp(param2set2zero(iii),'sigma_z')
        model_sol_new.param(indic_param) = -Inf;
        latexTable = [latexTable 'no cyclical growth &' ];
    end

    % avoid recalculation of endogenous parameters
    endo_flag = 0;
    model_sol_new = make_model_sol(model_sol_new,endo_flag);

    % Compute term premiums:
    [uncond_r_new,uncond_rn_new,uncond_TPr_new,uncond_TPrn_new,uncond_IRP_new,...
        all_stdv_r_new,all_stdv_rn_new,all_stdv_TPr_new,all_stdv_TPrn_new,all_stdv_IRP_new] =...
        compute_uncond_yds_TP(model_sol_new,H);

    % Nominal:
    for i = 1:nb_maturity
        latexTable = [latexTable ' & ' sprintf(Format,...
            all_stdv_TPrn_new(frequency * maturities_in_year(i)))];
    end
    latexTable = [latexTable ' & '];
    % Real:
    for i = 1:nb_maturity
        latexTable = [latexTable ' & ' sprintf(Format,...
            all_stdv_TPr_new(frequency * maturities_in_year(i)))];
    end
    latexTable = [latexTable ' & '];
    % Inflation:
    for i = 1:nb_maturity
        latexTable = [latexTable ' & ' sprintf(Format,...
            all_stdv_IRP_new(frequency * maturities_in_year(i)))];
    end

    latexTable = [latexTable ' \\ '];

end

latexTable = [latexTable '\hline \end{tabular}'];

if indic_save_output == 1
    % Save LaTeX table to a file
    latexFileName = 'Tables/table_decompTP.tex';
    fid = fopen(latexFileName, 'w');
    fprintf(fid, '%s', latexTable);
    fclose(fid);
    disp(['LaTeX table of descriptive statistics saved as ' latexFileName]);
end



