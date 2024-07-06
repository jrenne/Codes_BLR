function[model_sol] = compute_sdf(model)
% This function solves for the s.d.f. and derives the risk neutral dynamics
% of X_t.

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

lambda0 = ((1 - mu_gamma0) * delta * (Id - delta * Phi')^(-1) * Phi'...
    - mu_gamma0 * Id)* mu_c1;
lambda1 = - mu_gamma1 * mu_c1' * (delta * Phi * (Id - delta * Phi)^(-1) + Id);

eta0star = - log(delta) + mu_c0 + .5 * mu_c1' * Sigma * (Sigma') * mu_c1 +...
    lambda0' * Sigma * (Sigma') * mu_c1;
eta1star =  (lambda1 * Sigma * (Sigma') + Phi') * mu_c1 ;

muQ = Sigma * (Sigma') * lambda0;
PhiQ = Phi + Sigma * (Sigma') * lambda1';

%Sigma_1 = Sigma^(-1);
Sigma_1 = pinv(Sigma); % use pseudo-inverse for cases where dimension of X
% different from that of epsilon.

I = eye(n_Z);
Phi_ZXQ  = zeros(n_Z,n_X);
Phi_ZXXQ = zeros(n_Z,n_X^2);
for i = 1:n_Z
    J_i = kron(eye(n_eps),I(:,i)');
    Phi_ZXQ(i,:)  = model.Gamma0'*J_i'*Sigma_1*(PhiQ - Phi);
    Phi_ZXXQ(i,:) = reshape(model.Gamma1'*J_i'*Sigma_1*(PhiQ - Phi),1,n_X^2)';
end

model_sol          = model;
model_sol.eta0star = eta0star;
model_sol.eta1star = eta1star;
model_sol.lambda0  = lambda0;
model_sol.lambda1  = lambda1;
model_sol.muQ      = muQ;
model_sol.PhiQ     = PhiQ;
model_sol.muZQ     = kron((Sigma_1 * muQ)',eye(n_Z)) * model.Gamma0;
model_sol.Phi_ZXQ  = Phi_ZXQ + kron((Sigma_1 * muQ)',eye(n_Z)) * model.Gamma1;
model_sol.Phi_ZXXQ = Phi_ZXXQ;


