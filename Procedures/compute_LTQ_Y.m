function[PsiQ] = compute_LTQ_Y(model_sol,u)
% Computes the risk-neutral Laplace transform of Y

Sigma    = model_sol.Sigma;
Phi_Z    = model_sol.Phi_Z;
muQ      = model_sol.muQ;
PhiQ     = model_sol.PhiQ;
muZQ     = model_sol.muZQ;
Phi_ZXQ  = model_sol.Phi_ZXQ;
Phi_ZXXQ = model_sol.Phi_ZXXQ;
Gamma0   = model_sol.Gamma0;
Gamma1   = model_sol.Gamma1;

n_X   = size(model_sol.PhiQ,1);
n_Z   = size(model_sol.Phi_Z,1);
n_eps = size(model_sol.Sigma,2);

u_X  = u(1:n_X);
u_Z  = u((n_X+1):(n_X+n_Z));
u_XX = u((n_X+n_Z+1):end); % u_XX is of dimension (n_X * (n_X + 1)/2)

Mnx = make_Mnx(n_X);
Knx = make_Knx(n_X);

v0 = Sigma' * u_X + ...
    2 * kron(Sigma',muQ') * Mnx * u_XX + ...
    kron(eye(n_eps),u_Z') * Gamma0;
v1 = kron(eye(n_eps),u_Z') * Gamma1 + ...
    2 * Sigma' * reshape(Mnx * u_XX,n_X,n_X) * PhiQ;
V  = reshape(u_XX' * Mnx' * kron(Sigma,Sigma),n_eps,n_eps);

aux = (eye(n_eps) - 2 * V)^(-1);

Psi0   = u_X' * muQ + ...
    u_XX' * Mnx' * reshape(muQ * muQ',n_X^2,1) + ...
    u_Z' * muZQ + ...
    1/2 * log(det(aux)) + ...
    1/2 * v0' * aux * v0;
Psi_X  = PhiQ' * u_X + ...
    2 * kron(PhiQ',muQ') * Mnx * u_XX + ...
    Phi_ZXQ' * u_Z + ...
    v1' * aux * v0;
Psi_Z  = Phi_Z' * u_Z;
Psi_XX = Knx * (Phi_ZXXQ' * u_Z + ...
    kron(PhiQ',PhiQ') * Mnx * u_XX + ...
    1/2 * reshape(v1' * aux * v1,n_X^2,1));

PsiQ = struct;
PsiQ.Psi0    = Psi0;
PsiQ.Psi_X   = Psi_X;
PsiQ.Psi_Z   = Psi_Z;
PsiQ.Psi_XX  = Psi_XX;
PsiQ.Psi_Y   = [Psi_X;Psi_Z;Psi_XX];

