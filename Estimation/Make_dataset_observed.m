% =========================================================================
% Prepare dataset of observed variables and determine 
% std dev of measurement errors
% =========================================================================

% Create dataset based on names_of_variables:
Data_StateSpace.dataset = [];
for i = 1:size(names_of_variables,2)
    indic_variable = strcmp(data_names,names_of_variables{i});
    Data_StateSpace.dataset = [Data_StateSpace.dataset data(:,indic_variable)];
end

% ************ Delete certain data:

% indic_date = (dates=='01-Sep-2008') + (dates=='01-Dec-2008');
% Data_StateSpace.dataset(find(indic_date),2) = NaN;
% indic_var1 = find(strcmp(names_of_variables,'TIPSY02'));
% indic_var2 = find(strcmp(names_of_variables,'REALR20'));
% Data_StateSpace.dataset(find(indic_date),indic_var1:indic_var2) = NaN;

% dataset_copy = Data_StateSpace.dataset;
% 
% Data_StateSpace.dataset = nan(size(Data_StateSpace.dataset));
% Data_StateSpace.dataset(:,1:3) = dataset_copy(:,1:3); %macro
% Data_StateSpace.dataset(:,4:8) = dataset_copy(:,4:8); %nominal
% %Data_StateSpace.dataset(:,13:16) = dataset_copy(:,13:16); %survey

% ************ Determine standard deviation of measurement errors:

T = length(Data_StateSpace.dataset);

Data_StateSpace.stdv_measur = 0.1*ones(T,size(Data_StateSpace.dataset,2));
indic_variable = find(strcmp(names_of_variables,'Output gap'));
Data_StateSpace.stdv_measur(:,indic_variable) = 0.2;
indic_variable = contains(names_of_variables,'REALR');
Data_StateSpace.stdv_measur(dates<'01-Jan-1999',indic_variable) = 0.2; %dates<'01-Jan-2004'
indic_variable = find(strcmp(names_of_variables,'RTR'));
Data_StateSpace.stdv_measur(dates<'01-Jul-1984',indic_variable) = 0.2;

indic_date = or(dates=='01-Sep-2008',dates=='01-Dec-2008');
indic_variable = contains(names_of_variables,{'CPI inflation','TIPSY','REALR'});
Data_StateSpace.stdv_measur(indic_date,indic_variable) = 1.0;

% Plot fit obtained with Kalman filter:
% plot_KF;
