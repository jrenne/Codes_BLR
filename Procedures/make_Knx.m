function[M] = make_Knx(n)

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


K = zeros(n*(n+1)/2,n^2);
for i = 1:n^2
  M(S(i),i) = 1;
end


