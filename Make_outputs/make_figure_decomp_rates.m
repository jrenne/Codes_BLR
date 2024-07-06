% =========================================================================
% Figure showing rates decomposition
% =========================================================================

maturities_in_year = [0.25;2;5;10;20];

H = max(frequency*maturities_in_year);

% Compute model-implied term premiums:
[A,B,A4r,B4r] = compute_AB(model_sol,H);
[An,Bn,Cn,Dn,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,H);

% Compute model-implied yields:
real_yields = 100 * frequency * (ones(T,1)*A4r + X * B4r);
nom_yields  = 100 * frequency * (ones(T,1)*A4rn + X * B4rn + Z * C4rn + XX * D4rn);

% Compute specification of rates under Expectation Hypothesis (EH):
[E,V,A1,B1] = compute_EV(model_sol,1); % needed to compute expectations
[APr,BPr,APrn,BPrn] = compute_AB_real_nom_P(model_sol,H,...
    eta0r,eta1r,eta0rn,eta1rn,A1,B1);
% Compute model-implied yields under EH:
real_yields_P = 100 * frequency * (ones(T,1)*APr + all_xi_tt * BPr);
nom_yields_P  = 100 * frequency * (ones(T,1)*APrn + all_xi_tt * BPrn);

%IRP
E_Infl = nom_yields_P - real_yields_P; % Expected Inflation
inflation_RP =  nom_yields - real_yields - E_Infl;

% Real term premium
real_TP = real_yields - real_yields_P;

% Break-Even Inflation Rate
BEIR = nom_yields - real_yields; 

figure;
for m = 1:length(maturities_in_year)
    % nominal yields
    subplot(3,length(maturities_in_year),m);
    hold on;
    bar(dates,[real_yields_P(:,maturities_in_year(m)*frequency) real_TP(:,maturities_in_year(m)*frequency) E_Infl(:,maturities_in_year(m)*frequency) inflation_RP(:,maturities_in_year(m)*frequency)],'stacked');
    plot(dates,nom_yields(:,maturities_in_year(m)*frequency),'LineWidth',2);
    plot(dates,Data_StateSpace.dataset(:,3+m),'LineWidth',1.5,'Color','black');
    hold off;
    if m == length(maturities_in_year)
        legend('expected real rate','real term premium','expected inflation','inflation risk premium','fitted','observed');
    end
    title([num2str(maturities_in_year(m)) 'y nominal yield']);
    % TIPS yields
    subplot(3,length(maturities_in_year),length(maturities_in_year)+m);
    hold on;
    bar(dates,[real_yields_P(:,maturities_in_year(m)*frequency) real_TP(:,maturities_in_year(m)*frequency)],'stacked');
    plot(dates,real_yields(:,maturities_in_year(m)*frequency),'LineWidth',2);
    if m ~= 1
        plot(dates,Data_StateSpace.dataset(:,7+m),'LineWidth',1.5,'Color','black');
    end
    if m == length(maturities_in_year)
        plot(dates,Data_StateSpace.dataset(:,7+m),'LineWidth',1.5,'Color','black');
        legend('expected real rate','real term premium','fitted','observed');
    end
    hold off;
    title([num2str(maturities_in_year(m)) 'y TIPS yield']);

    % BEIR
    subplot(3,length(maturities_in_year),2*length(maturities_in_year)+m);
    hold on;
    bar(dates,[E_Infl(:,maturities_in_year(m)*frequency) inflation_RP(:,maturities_in_year(m)*frequency)],'stacked');
    plot(dates,BEIR(:,maturities_in_year(m)*frequency),'LineWidth',2);
    if m ~=1
        plot(dates,Data_StateSpace.dataset(:,3+m)-Data_StateSpace.dataset(:,7+m),'LineWidth',1.5,'Color','black');
    end
    if m == length(maturities_in_year)
        plot(dates,Data_StateSpace.dataset(:,3+m)-Data_StateSpace.dataset(:,7+m),'LineWidth',1.5,'Color','black');
        legend('expected inflation','inflation risk premium','fitted','observed');
    end
    hold off;
    title([num2str(maturities_in_year(m)) 'y BEIR']);
end
