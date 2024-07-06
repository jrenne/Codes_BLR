% =========================================================================
% Create Latex table of model parameterization
% =========================================================================

T = length(Data_StateSpace.dataset);
K = length(model_sol.param);

alpha = 5/100;

list_param_in_table = model_sol.names_param;

nb_param = size(list_param_in_table,1);

param_desc = cell(nb_param,1);
param_values = zeros(nb_param,1);
param_stdev  = zeros(nb_param,1);
param_names  = cell(nb_param,1);

param_desc{1}  = 'AR trend growth';
param_desc{2}  = 'AR consumption gap';
param_desc{3}  = 'AR fast risk aversion';
param_desc{4}  = 'AR slow risk aversion';
param_desc{5}  = 'AR $corr(\pi,c)$ factor';
param_desc{6}  = 'SD trend growth';
param_desc{7}  = 'SD output gap';
param_desc{8}  = 'SD fast risk aversion';
param_desc{9}  = 'SD slow risk aversion';
param_desc{10} = 'SD $corr(\pi,c)$ factor';
param_desc{11} = 'AR trend inflation';
param_desc{12} = 'SD trend inflation';
param_desc{13} = 'AR cyclical inflation';
param_desc{14} = 'SD cyclical inflation';
param_desc{15} = 'mean inflation';
param_desc{16} = 'mean consumption growth';
param_desc{17} = 'mean risk aversion';
param_desc{18} = 'mean $corr(\pi,c)$ factor';
param_desc{19} = 'time discount factor';
param_desc{20} = 'hysteresis effect';

for i = 1:nb_param
    indic_variable  = strcmp(model_sol.names_param,list_param_in_table{i});
    param_values(i) = model_sol.param_transf(indic_variable);
    param_stdev(i)  = param_transf_stdev(indic_variable);
    param_names(i)  = model_sol.names_param_Latex(indic_variable);
end

latexTable = sprintf('\\begin{tabular}{llrrrr} \\hline ');
latexTable = [latexTable 'Description & Parameter & Point estimate & Std. Dev. & \multicolumn{2}{c}{Confidence range} \\ \hline '];

for i = 1:nb_param
    latexTable = [latexTable param_desc{i}];
    latexTable = [latexTable ' & ' param_names{i}];
    latexTable = [latexTable ' & $' sprintf('%0.5f',param_values(i)) '$'];
    latexTable = [latexTable ' & $' sprintf('%0.5f',param_stdev(i))  '$'];
    latexTable = [latexTable ' & $' sprintf('%0.5f',param_values(i) - tinv(1-alpha,T-K)*param_stdev(i))  '$'];
    latexTable = [latexTable ' & $' sprintf('%0.5f',param_values(i) + tinv(1-alpha,T-K)*param_stdev(i))  '$'];
    latexTable = [latexTable ' \\ '];
end

latexTable = [latexTable ' \hline '];
latexTable = [latexTable ' \end{tabular}'];

% param_per_line = 5;
% nb_of_lines = ceil(nb_param/param_per_line);
% 
% latexTable = sprintf('\\begin{tabular}{c');
% for i = 2:param_per_line
%     latexTable = [latexTable 'c'];
% end
% 
% latexTable = [latexTable '} \\ \hline'];
% 
% count = 0;
% for row = 1:nb_of_lines
%     % name of variables
%     count_before_line = count;
% 
%     % Line with names of variables:
%     count_this_line = 0;
%     for j = 1:param_per_line
%         count_this_line = count_this_line + 1;
%         count = count + 1;
%         if count == nb_param
%             latexTable = [latexTable param_names{count} '\\'];
%             break;
%         else
%             if count_this_line < param_per_line
%                 latexTable = [latexTable param_names{count} '&'];
%             else
%                 latexTable = [latexTable param_names{count} '\\'];
%             end
%         end
%     end
% 
%     % Line with values of variables:
%     count = count_before_line;
%     count_this_line = 0;
%     for j = 1:param_per_line
%         count_this_line = count_this_line + 1;
%         count = count + 1;
%         if count == nb_param
%                 latexTable = [latexTable '$' sprintf('%0.3g',param_values(count)) '$\\ \hline'];
%                 break
%         else
%             if count_this_line < param_per_line
%                 latexTable = [latexTable '$' sprintf('%0.3g',param_values(count)) '$&'];
%             else
%                 latexTable = [latexTable '$' sprintf('%0.3g',param_values(count)) '$\\ \hline'];
%             end
%         end
%     end
% end
% 
% %latexTable = [latexTable '\\ \hline'];
% 
% latexTable = [latexTable ' \end{tabular}'];

if indic_save_output == 1
    % Save LaTeX table to a file
    latexFileName = 'Tables/table_param.tex';
    fid = fopen(latexFileName, 'w');
    fprintf(fid, '%s', latexTable);
    fclose(fid);
    disp(['LaTeX table of descriptive statistics saved as ' latexFileName]);
end

