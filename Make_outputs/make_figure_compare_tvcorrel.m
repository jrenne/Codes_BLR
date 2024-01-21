% Compare model-based pi-cons correlation

myLW=2;
disp('----------------------');
disp('analyse model implied pi cons correlation with that of time series data');
disp('----------------------');
WinTvCorr=20;
disp(['tv correlation window (quarters): ',num2str(WinTvCorr)]);
GrVarsInds=[1 2 25 38];  % indices of growth variables
InflVarsInds=[4 5 6 7 8]; % indices of inflation variables
nbGr=length(GrVarsInds);
nbInfl=length(InflVarsInds);

disp('Candidate growth variables:');
disp(data_names(GrVarsInds)');
disp('Candidate inflation variables:');
disp(data_names(InflVarsInds)');
GrVars=data(:,GrVarsInds);
GrVarsNames=data_names(GrVarsInds)';
InflVars=data(:,InflVarsInds);
InflVarsNames=data_names(InflVarsInds)';
[COEFF, SCORE] = pca(GrVars); PCGrVars=SCORE(:,1);
[COEFF, SCORE] = pca(InflVars); PCInflVars=SCORE(:,1);
GrVars=[GrVars PCGrVars];
InflVars=[InflVars PCInflVars];

TT=length(dates);
ctr=0;
TvCorrAll=NaN(TT,(nbGr+1)*(nbInfl+1));
IndCatch=NaN((nbGr+1)*(nbInfl+1),3);
for gg=1:nbGr+1;
    for ii=1:nbInfl+1;
        ctr=ctr+1;
        IndCatch(ctr,:)=[ctr gg ii];
        GrLoc=GrVars(:,gg);
        InflLoc=InflVars(:,ii);
        tvcorrLoc=NaN(TT,1);
        for tt=WinTvCorr:TT;
            tvcorrLoc(tt)=corr(GrLoc(tt-WinTvCorr+1:tt),InflLoc(tt-WinTvCorr+1:tt));
        end;
        TvCorrAll(:,ctr)=tvcorrLoc;
    end;
end;
%[HPTrTvCorrAll,Cyclical] = hpfilter(TvCorrAll(WinTvCorr:end,:),1600);
%HPTrTvCorrAll=[NaN(WinTvCorr-1,(nbGr+1)*(nbInfl+1)); HPTrTvCorrAll];

[HPTrTvCorrAll,Cyclical] = hpfilter(TvCorrAll(WinTvCorr:end,1:18),1600);
HPTrTvCorrAll=[NaN(WinTvCorr-1,18); HPTrTvCorrAll];

% find max correlation 
CorCatch=corr(HPTrTvCorrAll(WinTvCorr:end,:),condCorrel(WinTvCorr:end,2));
[SortCorCatch,SortCorCatchInd]=sort(CorCatch);
maxInd=SortCorCatchInd(end);
HPTrTvCorrBest=HPTrTvCorrAll(:,maxInd);

disp('Highest correlation between HP trend tv correlation and model correlation is obtained for the following pair');
disp(['Growth: ',GrVarsNames{IndCatch(maxInd,2)}]);
disp(['Inflation: ',InflVarsNames{IndCatch(maxInd,3)}]);

figure(401);
plot(dates,TvCorrAll);

figure(402);
plot(dates,HPTrTvCorrAll);

figure(403);
hold on;
title('\pi vs cons correlation');
plot(dates,condCorrel(:,2),'LineWidth',myLW);
plot(dates,HPTrTvCorrBest,'LineWidth',myLW);
legend({'Model','emp time series'}); 
