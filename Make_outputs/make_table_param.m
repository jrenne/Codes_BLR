% =========================================================================
% Create Latex table of model parameterization
% =========================================================================

param_per_line = 8;

model_sol = make_model_sol(model_sol);
param2print = array2table(round(model_sol.param_transf,4)',...
    'RowNames',model_sol.names_param);
disp(param2print);

list_param_in_table = {'rho_g';'rho_z';'rho_w';'rho_k';...
    'sigma_g';'sigma_z';'sigma_w';'sigma_k';...
    'mu_gamma1_g';'mu_gamma1_z';'rho_pi';'rho_pi_star';...
    'sigma_pi_star';'sigma_pi_z';'mu_pi';'mu_c';...
    'mu_gamma0';'delta';'rho_gz'};

nb_param = size(list_param_in_table,1);

param_values = zeros(nb_param,1);
param_names  = cell(nb_param,1);

for i = 1:nb_param
    indic_variable = find(strcmp(model_sol.names_param,list_param_in_table{i}));
    param_values(i) = model_sol.param_transf(indic_variable);
    param_names(i)  = model_sol.names_param_Latex(indic_variable);
end


nb_of_lines = ceil(nb_param/param_per_line);

latexTable = sprintf('\\begin{tabular}{c');
for i = 2:param_per_line
    latexTable = [latexTable 'c'];
end

latexTable = [latexTable '} \\ \hline'];

count = 0;
for row = 1:nb_of_lines
    % name of variables
    count_before_line = count;

    % Line with names of variables:
    count_this_line = 0;
    for j = 1:param_per_line;
        count_this_line = count_this_line + 1;
        count = count + 1;
        if count == nb_param
            latexTable = [latexTable param_names{count} '\\'];
            break;
        else
            if count_this_line < param_per_line
                latexTable = [latexTable param_names{count} '&'];
            else
                latexTable = [latexTable param_names{count} '\\'];
            end
        end
    end

    % Line with values of variables:
    count = count_before_line;
    count_this_line = 0;
    for j = 1:param_per_line;
        count_this_line = count_this_line + 1;
        count = count + 1;
        if count == nb_param
                latexTable = [latexTable '$' sprintf('%0.3g',param_values(count)) '$\\ \hline'];
                break
        else
            if count_this_line < param_per_line
                latexTable = [latexTable '$' sprintf('%0.3g',param_values(count)) '$&'];
            else
                latexTable = [latexTable '$' sprintf('%0.3g',param_values(count)) '$\\ \hline'];
            end
        end
    end
end

%latexTable = [latexTable '\\ \hline'];

latexTable = [latexTable ' \end{tabular}'];


% Save LaTeX table to a file
latexFileName = 'Tables/table_param.tex';
fid = fopen(latexFileName, 'w');
fprintf(fid, '%s', latexTable);
fclose(fid);

disp(['LaTeX table of descriptive statistics saved as ' latexFileName]);

