function vechX = vech(X)

n = size(X,1);
aux = reshape(1:n^2,n,n);
indicators = nonzeros(tril(aux));

vechX = X(indicators);


