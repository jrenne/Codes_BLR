function[A,B,A4r,B4r] = compute_AB(model_sol,H)
% This functions computes the term structure of real rates.

eta0star = model_sol.eta0star;
eta1star = model_sol.eta1star;
muQ      = model_sol.muQ;
PhiQ     = model_sol.PhiQ;
Sigma    = model_sol.Sigma;

n_X = size(PhiQ,1);

A = zeros(1,H);
B = zeros(n_X,H);
A4r = zeros(1,H);
B4r = zeros(n_X,H);

a = -eta0star;
b = -eta1star;

A(:,1) = a;
B(:,1) = b;
A4r(:,1) = - a;
B4r(:,1) = - b;

for h = 2:H
    a = - eta0star + a + b' * muQ + .5 * b' * Sigma * (Sigma') * b;
    b = - eta1star + PhiQ' * b ;
    A(:,h) = a;
    B(:,h) = b;
    A4r(:,h) = -a/h;
    B4r(:,h) = -b/h;
end

