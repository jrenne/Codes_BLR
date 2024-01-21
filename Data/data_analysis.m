%clear
%close all
%clc

global frequency;

% Choose which data file to load
frequency = 4; % 12: monthly, 4: quarterly

if frequency == 12
    [X,TXT,~] = xlsread('Data/US_data_matlab.xlsx','monthly');
    X(:,2:3) = X(:,2:3)./repmat(X(:,4),1,2); % calculate per-capita consumption
    X = X(1:end-24,:); % stop in Dec 2019
elseif frequency == 4
    [X,TXT,~] = xlsread('Data/US_data_matlab.xlsx','quarterly');
    X(:,2:3) = X(:,2:3)./repmat(X(:,4),1,2); % calculate per-capita consumption
    X(:,26:27) = X(:,26:27)./repmat(X(:,4),1,2); % calculate per-capita GDP
    X = X(1:end-8,:); % stop in Q4 2019
end

data_names = TXT(1,2:end); % save variable names
dates = x2mdate(X(2:end,1),0,'datetime'); % save date column
data = (log(X(2:end,2:9)) - log(X(1:end-1,2:9)))*100; % calculate annualised period-on-period growth rates
data(:,9:24) = X(2:end,10:25); % yields
data(:,25:26)= (log(X(2:end,26:27)) - log(X(1:end-1,26:27)))*100; % real GDP (+ potential)

% SPF data:
data(:,27:33) = X(2:end,28:34);

[num_obs, num_var] = size(data); % determine sample size (observations x variables)

% Calculate slopes:
data(:,end+1) = data(:,15) - data(:,11); % nominal 10y - 2y 
data(:,end+1) = data(:,22) - data(:,18); % real 10y - 2y
data(:,end+1) = data(:,end-1) - data(:,end); % BEIR 10y - 2y 
data_names(end+1) = {'nominal slope'};
data_names(end+1) = {'real slope'};
data_names(end+1) = {'BEIR slope'};
num_var = num_var + 3;

% Add output gap:
data(:,end+1) = log(X(2:end,26)./X(2:end,27));
data_names(end+1) = {'output gap'};
num_var = num_var + 1;

% Add 10-year term premium:
data(:,end+1) = data(:,15) - data(:,30);
data_names(end+1) = {'10-year nom. term premium'};
num_var = num_var + 1;

% Plot autocorrelation functions
data_acf = nan(21,num_var);
figure('Name','Autocorrelations','WindowState','maximized');
for v = 1:num_var
    subplot(7,7,v);
    [data_acf(:,v),~,~,~] = autocorr(data(:,v));
    title(data_names{v});
end

% Plot all data
figure('Name','Data','WindowState','maximized');
for v = 1:num_var
    subplot(7,7,v);
    % detect if sparse data:
    sparse_data = (sum((~isnan(data(2:end,v))).*(~isnan(data(1:(end-1),v))))==0);
    if sparse_data
        plot(dates,data(:,v),'o');
    else
        plot(dates,data(:,v));
    end
    title(data_names{v});
    xlim([dates(1) dates(end)]);
end

% % Plot average yield curves
% figure('Name','Term structures','WindowState','maximized');
% hold on;
% plot([3/12,1,2,3,4,5,10,15,20],mean(data(~isnan(data(:,18)),9:17),'omitnan'));
% plot([2,3,4,5,10,15,20],mean(data(:,18:24),'omitnan'));
% legend('nominal average since 1999','real average since 1999');
% xlabel('years to maturity');
% hold off;

% Calculate and display data moments
data_mean = mean(data,'omitnan');
data_std = std(data,'omitnan');
data_corr = corr(data,'rows','pairwise'); %pairwise correlations ignoring NAN-values (possibly not positive semi-definite!)

% disp(array2table(round([sum(~isnan(data))', data_mean', data_std', data_acf(2:4,:)'],4),'RowNames',data_names,'VariableNames',{'# obs','mean', 'std', 'ac(1)', 'ac(2)', 'ac(3)'}));
% disp('Correlation matrix:');
% disp(array2table(round(data_corr,2),'RowNames',string(1:num_var),'VariableNames',string(1:num_var)));


%% Tests with dates:
%DateStrings = {'2014-05-26';'2014-08-03'};
%t = datetime(DateStrings,'InputFormat','yyyy-MM-dd')

% Sample recession periods (replace with actual recession periods)
recessionDates = [
    datetime('01-Jan-1970'), datetime('30-Nov-1970');
    datetime('01-Dec-1973'), datetime('31-Mar-1975');
    datetime('01-Jan-1980'), datetime('31-Jul-1980');
    datetime('01-Jul-1981'), datetime('30-Nov-1982');
    datetime('01-Jul-1990'), datetime('31-Mar-1991');
    datetime('01-Mar-2001'), datetime('30-Nov-2001');
    datetime('01-Dec-2007'), datetime('30-Jun-2009');
];
recessionStarts = recessionDates(:,1);
recessionEnds   = recessionDates(:,2);



% Create dataset where nominal yields are NaN when it is the case for
% {'TIPSY10'}:
indic_TIPSnan = find(isnan(data(:,find(strcmp(data_names,{'TIPSY10'})))));
indic_fst_nom_yield = find(strcmp(data_names,{'YIELD3M'}));
indic_lst_nom_yield = find(strcmp(data_names,{'YIELD20'}));
data_bis = data;
data_bis(indic_TIPSnan,indic_fst_nom_yield:indic_lst_nom_yield) = nan;
indic_nom_slope = find(strcmp(data_names,{'nominal slope'}));
data_bis(indic_TIPSnan,indic_nom_slope) = nan;
indic_inflation = find(strcmp(data_names,{'CPI all '}));
data_bis(indic_TIPSnan,indic_inflation) = nan;
