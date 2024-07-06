% =========================================================================
% Create Latex table of summary statistics
% =========================================================================

% Convert the data matrix to a table
dataTable = array2table(Data_StateSpace.dataset,...
    'VariableNames', Data_StateSpace.names_of_variables);

%Format = '%.4f'; % format of numbers

% Compute descriptive statistics
meanValues   = nanmean(Data_StateSpace.dataset);
stdDevValues = nanstd(Data_StateSpace.dataset);
minValues    = min(Data_StateSpace.dataset);
maxValues    = max(Data_StateSpace.dataset);

% detect when variable become available:
nb_var = min(size(Data_StateSpace.dataset));
dates_start = strings(nb_var,1);
for i = 1:nb_var
    dates_start(i) =  datestr(dates(find(~isnan(Data_StateSpace.dataset(:,i)), 1)));
end

% define variable type
variable_type = cell(1,length(Data_StateSpace.names_of_variables));
variable_type{1}  = 'Macro';
variable_type{4}  = 'Nominal';
variable_type{5}  = 'yield';
variable_type{9}  = 'Real';
variable_type{10} = 'yield';
variable_type{13} = 'Survey';

% define model counterpart
variable_counterpart = cell(1,length(Data_StateSpace.names_of_variables));
variable_counterpart{1}  = '$\Delta c_t$';
variable_counterpart{2}  = '$\pi_t$';
variable_counterpart{3}  = '$z_t$';
variable_counterpart{4}  = '$i_t$';
variable_counterpart{5}  = '$i_{t,8}$';
variable_counterpart{6}  = '$i_{t,20}$';
variable_counterpart{7}  = '$i_{t,40}$';
variable_counterpart{8}  = '$i_{t,80}$';
variable_counterpart{9}  = '$r_{t,8}$';
variable_counterpart{10} = '$r_{t,20}$';
variable_counterpart{11} = '$r_{t,40}$';
variable_counterpart{12} = '$r_{t,80}$';
variable_counterpart{13} = '$\mathbb{E}_t i_{t+40}$';
variable_counterpart{14} = '$\mu_\pi^{PCE} + 4\mathbb{E}_t \pi_{t+40}^*$';
variable_counterpart{15} = '$\frac{1}{40}\sum_{h=1}^{40}\mathbb{E}_t i_{t+h}$';
variable_counterpart{16} = '$\frac{4}{40}\sum_{h=1}^{40}\mathbb{E}_t \pi_{t+h}$';

% Create a cell array for the table content
tableContent = [Data_StateSpace.names_of_variables; ...
                cellfun(@num2str, num2cell(meanValues), 'UniformOutput', false); ...
                cellfun(@num2str, num2cell(stdDevValues), 'UniformOutput', false); ...
                cellfun(@num2str, num2cell(minValues), 'UniformOutput', false); ...
                cellfun(@num2str, num2cell(maxValues), 'UniformOutput', false)];

% Shorten long name
tableContent{1,1} = 'Consumption';

% Create LaTeX table
latexTable = sprintf('\\begin{tabular}{llrrrrrc} \\hline ');
latexTable = [latexTable 'Type & Variable & Mean & S.D. & Min & Max & First date & Model $\times 100$ \\ \hline '];

for row = 1:size(tableContent, 2)
    latexTable = [latexTable variable_type{1, row}];
    latexTable = [latexTable ' & ' tableContent{1, row}];
    for col = 2:size(tableContent,1)
    latexTable = [latexTable ' & $' sprintf('%.2f',str2num(tableContent{col, row})) '$'];
    end
    latexTable = [latexTable ' & ' extractAfter(dates_start(row),"01-")];
    latexTable = [latexTable ' & ' variable_counterpart{1, row}];
    latexTable = [latexTable ' \\ '];
end

latexTable = [latexTable '\hline \end{tabular}'];

if indic_save_output == 1
    % Save LaTeX table to a file
    latexFileName = 'Tables/table_desc_stat.tex';
    fid = fopen(latexFileName, 'w');
    fprintf(fid, '%s', latexTable);
    fclose(fid);
    disp(['LaTeX table of descriptive statistics saved as ' latexFileName]);
end
