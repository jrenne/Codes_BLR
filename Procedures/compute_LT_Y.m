function[Psi] = compute_LT_Y(model,u)
% Computes the physical Laplace transform of Y

Sigma = model.Sigma;
Phi   = model.Phi;
Phi_Z = model.Phi_Z;

n_X   = size(model.Phi,1);
n_Z   = size(model.Phi_Z,1);
n_eps = size(model.Sigma,2);

u_X  = u(1:n_X);
u_Z  = u((n_X+1):(n_X+n_Z));
u_XX = u((n_X+n_Z+1):end); % u_XX is of dimension (n_X * (n_X + 1)/2)

Mnx = make_Mnx(n_X);
Knx = make_Knx(n_X);

v0 = Sigma' * u_X + kron(eye(n_eps),u_Z') * model.Gamma0;
v1 = kron(eye(n_eps),u_Z') * model.Gamma1 + ...
    2 * Sigma' * reshape(Mnx * u_XX,n_X,n_X) * Phi;
V  = reshape(u_XX' * Mnx' * kron(Sigma,Sigma),n_eps,n_eps);

aux = (eye(n_eps) - 2 * V)^(-1);

Psi0   = 1/2 * log(det(aux)) + 1/2 * v0' * aux * v0;
Psi_X  = Phi' * u_X + v1' * aux * v0;
Psi_Z  = Phi_Z' * u_Z;
Psi_XX = Knx * (kron(Phi',Phi') * (Mnx * u_XX) + ...
    1/2 * reshape(v1' * aux * v1,n_X^2,1));

Psi = struct;
Psi.Psi0    = Psi0;
Psi.Psi_X   = Psi_X;
Psi.Psi_Z   = Psi_Z;
Psi.Psi_XX  = Psi_XX;
Psi.Psi_Y   = [Psi_X;Psi_Z;Psi_XX];

