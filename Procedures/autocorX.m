function[acn] = autocorX(y,n)

%[aux,~,~,~] = autocorr(y,n);
%acn = aux(n);

X = y((n+1):end);
Y = y(1:(end-n));

aux = [X Y];
aux = aux(~isnan(X+Y),:);
aux = corr(aux);
acn = aux(1,2);
