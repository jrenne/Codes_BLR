function[Delta_c,gamma,inflation,X,Z,XX,Y] = simul_model(model,T,x0,z0)
% This function simulates the model on T periods, starting from x0.

n_X   = size(model.Phi,1);
n_Z   = size(model.Phi_Z,1);
n_eps = size(model.Sigma,2);

epsilon = normrnd(0,1,n_eps,T);

x = x0;

X = x0';
XX = vech(x0 * x0')';

z = z0;
Z = z0';

for t = 2:T
    
    x = model.Phi * x + model.Sigma * epsilon(:,t);
    
    SigmaZ = reshape(model.Gamma0 + model.Gamma1 * x,n_Z,n_eps);
    
    %disp(SigmaZ);
    
    z      = model.Phi_Z * z + SigmaZ * epsilon(:,t);
    
    X  = [X;x'];
    Z  = [Z;z'];
    XX = [XX;vech(x * x')'];
end

Delta_c   = model.mu_c0 + X * model.mu_c1;
gamma     = model.mu_gamma0 + X * model.mu_gamma1;
inflation = model.mu_pi0 + ...
    Z * model.mu_piZ + ...
    X * model.mu_piX + ...
    XX * model.mu_piXX;

Y = [X Z XX];
