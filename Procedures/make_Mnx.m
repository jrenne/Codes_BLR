function[M] = make_Mnx(n)

S = zeros(n,n);
aux = 1:n*(n+1)/2;
count = 0;
for j = 1:n
    nb = n - j + 1;
    S(j:n,j) = aux((count+1):(count + nb));
    count = count + nb;
end
count = 0;
for j = 1:n
    nb = n - j + 1;
    S(j,(j+1):n) = aux((count+2):(count + nb));
    count = count + nb;
end


M = zeros(n^2,n*(n+1)/2);
for i = 1:n^2
  M(i,S(i)) = .5;
end

diagS = diag(S);
G = reshape(1:n^2,n,n);
diagG = diag(G);
for i = 1:n
  M(diagG(i),diagS(i)) = 1;
end




