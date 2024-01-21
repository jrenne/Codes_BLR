function[X,Z,XX,Y] = simul_model_multi(model,T,x0,z0,nb_replic)
% This function simulates the model on T periods, starting from x0.

n_X   = size(model.Phi,1);
n_Z   = size(model.Phi_Z,1);
n_eps = size(model.Sigma,2);

epsilon = normrnd(0,1,n_eps*T*nb_replic,1);
epsilon = reshape(epsilon,n_eps,nb_replic,T);

x = x0 * ones(1,nb_replic);
X = zeros(n_X,nb_replic,T);
X(:,:,1) = x;

z = z0 * ones(1,nb_replic);
Z = zeros(n_Z,nb_replic,T);
Z(:,:,1) = z;

XX = zeros(n_X^2,nb_replic,T);
XX(:,:,1) = kron(x0,x0) * ones(1,nb_replic);

for t = 2:T
    
    x = model.Phi * x + model.Sigma * epsilon(:,:,t);
    
    SigmaZ = model.Gamma0 + model.Gamma1 * x;
    
    SigmaZeps = SigmaZ .* kron(epsilon(:,:,t),ones(n_Z,1));
    SigmaZeps = kron(ones(1,n_eps),eye(n_Z)) * SigmaZeps;
    
    z      = model.Phi_Z * z + SigmaZeps;
    
    X(:,:,t)  = x;
    Z(:,:,t)  = z;
    XX(:,:,t) = kron(x,ones(n_X,1)) .* kron(ones(n_X,1),x);
    
end

% Delta_c   = model.mu_c0 + X * model.mu_c1;
% gamma     = model.mu_gamma0 + X * model.mu_gamma1;
% inflation = model.mu_pi0 + ...
%     Z * model.mu_piZ + ...
%     X * model.mu_piX + ...
%     XX * model.mu_piXX;

Y = [X;Z;XX];
