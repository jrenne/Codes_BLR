function[u0,u1] = compute_u(model,X0,X1)
% This function solves for the utility function.
% It is essentially used to check the formulas.

delta     = model.delta;
mu_gamma0 = model.mu_gamma0;
mu_gamma1 = model.mu_gamma1;
mu_c0     = model.mu_c0;
mu_c1     = model.mu_c1;
Sigma     = model.Sigma;
Phi       = model.Phi;

n_X   = size(model.Phi,1);
n_Z   = size(model.Phi_Z,1);
n_eps = size(model.Sigma,2);

Id = eye(n_X);

mu_u1  = delta * (Id - delta * Phi')^(-1) * (Phi') * mu_c1;

mu_u0t   = 0;

gammat = mu_gamma0 + mu_gamma1' * X0;
aux = (ones(n_X,1) * (1-gammat)) .* (((Id - delta * Phi')^(-1) * Phi' + Id) * mu_c1);
disp(size(aux));

aux = (aux') * Sigma;
aux = kron(aux,ones(1,n_eps)) .* kron(ones(1,n_eps),aux);
aux = aux * ones(n_eps^2,1);

mu_u0tp1 = - mu_c0 + 1/delta * ...
(mu_u0t - .5 * delta * (aux')./(1-gammat));

u0 =                         mu_u0t   + (mu_u1') * X0;
u1 = (mu_c0 + mu_c1' * X1) + mu_u0tp1 + (mu_u1') * X1;
