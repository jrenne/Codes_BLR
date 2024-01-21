
% Define parameters to be optimized on:
FILTER = 0 * model_sol.param + 1;
% Create vector of parameters:
sub_parameters =  model_sol.param(FILTER==1);

% Compute moments:
indic_disp = 0;
M = print_moments(sub_parameters,model_sol,indic_disp);


rowNames = M.Properties.RowNames;
colNames = M.Properties.VariableNames(1:4);

% Create LaTeX table
latexTable = sprintf('\\begin{tabular}{lrrrr} \\hline ');
for i = 1:size(colNames,2)
    latexTable = [latexTable ' & ' colNames{i}];
end
latexTable = [latexTable ' \\ \hline '];

for row = 1:size(rowNames, 1)
    latexTable = [latexTable rowNames{row}];
    for col = 1:size(colNames,2)
        if col < size(colNames,2) % not the last entry
            if abs(M{row, col}) < 1
                Format = '%.4f';
            elseif abs(M{row, col}) >= 100
                Format = '%.0f';
            else
                Format = '%.2f';
            end
        else % last entry of column (i.e. weight)
            if abs(M{row, col}) < 1
                Format = '%.2f';
            else
                Format = '%.0f';
            end
        end
        latexTable = [latexTable ' & $' sprintf(Format,M{row, col}) '$'];
    end
    latexTable = [latexTable ' \\ '];
end

latexTable = [latexTable '\hline \end{tabular}'];

% Save LaTeX table to a file
latexFileName = 'Tables/table_moments.tex';
fid = fopen(latexFileName,'w');
fprintf(fid, '%s', latexTable);
fclose(fid);

disp(['LaTeX table of descriptive statistics saved as ' latexFileName]);
