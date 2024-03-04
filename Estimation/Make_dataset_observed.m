% =========================================================================
% Prepare dataset of observed variables and determine 
% std dev of measurement errors
% =========================================================================

% Create dataset based on names_of_variables:
Data_StateSpace.dataset = [];
for i = 1:size(names_of_variables,2)
    indic_variable = find(strcmp(data_names,names_of_variables{i}));
    Data_StateSpace.dataset = [Data_StateSpace.dataset data(:,indic_variable)];
end

%names_of_variables2interpolate = {'RGDP10','BOND10','BILL10'};
names_of_variables2interpolate = {'BOND10','BILL10'};
for i = 1:size(names_of_variables2interpolate,2)
    indic_variable = find(strcmp(names_of_variables,names_of_variables2interpolate{i}));
    x = Data_StateSpace.dataset(:,indic_variable);
    nanx = isnan(x);
    t    = 1:numel(x);
    x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
    Data_StateSpace.dataset(:,indic_variable) = x;
end


% Determine measurement errors:
Data_StateSpace.stdv_measur = sqrt(var(Data_StateSpace.dataset,"omitnan"))/4;
%Data_StateSpace.stdv_measur(1:2) = 0.05; % for consumption growth and inflation
Data_StateSpace.stdv_measur(1:2) = 0.2; % for consumption growth and inflation

% Low std dev for fit of BILL10:
indic_variable = find(strcmp(names_of_variables,'BILL10'));
%Data_StateSpace.stdv_measur(indic_variable) = .1 ;
Data_StateSpace.stdv_measur(indic_variable) = .2 ;

% Low std dev for fit of nominal yields:
Data_StateSpace.stdv_measur(4:6) = .2 ;
%Data_StateSpace.stdv_measur(4:6) = .05 ;

% Low std dev for fit of real yields:
Data_StateSpace.stdv_measur(7:8) = .2 ;

% Plot fit obtained with Kalman filter:
% plot_KF;
