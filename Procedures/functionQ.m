function[Q] = functionQ(StateSpace,xi)
% This function computes the conditional variance matrix of the latent
% variable

n_xi = size(StateSpace.F,1);

Theta0 = StateSpace.Theta0;
Theta1 = StateSpace.Theta1;

Q = reshape(Theta0 + Theta1 * xi,n_xi,n_xi);



