% =========================================================================
% Create Latex table of summary statistics
% =========================================================================

% Convert the data matrix to a table
dataTable = array2table(Data_StateSpace.dataset,...
    'VariableNames', Data_StateSpace.names_of_variables);

%Format = '%.4f'; % format of numbers

% Compute descriptive statistics
meanValues = nanmean(Data_StateSpace.dataset);
stdDevValues = nanstd(Data_StateSpace.dataset);
minValues = min(Data_StateSpace.dataset);
maxValues = max(Data_StateSpace.dataset);

% detect when variable become available:
nb_var = min(size(Data_StateSpace.dataset));
dates_start = strings(nb_var,1);
for i = 1:nb_var
    dates_start(i) =  datestr(dates(find(~isnan(Data_StateSpace.dataset(:,i)), 1)));
end


% Create a cell array for the table content
tableContent = [Data_StateSpace.names_of_variables; ...
                cellfun(@num2str, num2cell(meanValues), 'UniformOutput', false); ...
                cellfun(@num2str, num2cell(stdDevValues), 'UniformOutput', false); ...
                cellfun(@num2str, num2cell(minValues), 'UniformOutput', false); ...
                cellfun(@num2str, num2cell(maxValues), 'UniformOutput', false)];

% Create LaTeX table
latexTable = sprintf('\\begin{tabular}{lrrrrr} \\hline ');
latexTable = [latexTable 'Variable & Mean & Std. Dev. & Min & Max & First date \\ \hline '];

for row = 1:size(tableContent, 2)
    latexTable = [latexTable tableContent{1, row}];
    for col = 2:size(tableContent,1)
    latexTable = [latexTable ' & $' sprintf('%.2f',str2num(tableContent{col, row})) '$'];
    end
    latexTable = [latexTable ' & ' extractAfter(dates_start(row),"01-")];
    latexTable = [latexTable ' \\ '];
end

latexTable = [latexTable '\hline \end{tabular}'];

% Save LaTeX table to a file
latexFileName = 'Tables/table_desc_stat.tex';
fid = fopen(latexFileName, 'w');
fprintf(fid, '%s', latexTable);
fclose(fid);

disp(['LaTeX table of descriptive statistics saved as ' latexFileName]);

