function[A,B,C,D,A4r,B4r,C4r,D4r] = compute_AB_nominal(model_sol,H)
% This functions computes the term structure of nominal rates.

eta0star = model_sol.eta0star;
eta1star = model_sol.eta1star;

mu_pi0  = model_sol.mu_pi0;
mu_piZ  = model_sol.mu_piZ;
mu_piX  = model_sol.mu_piX;
mu_piXX = model_sol.mu_piXX;

n_X   = size(model_sol.PhiQ,1);
n_Z   = size(model_sol.Phi_Z,1);

A   = zeros(1,H);
B   = zeros(n_X,H);
C   = zeros(n_Z,H);
D   = zeros(n_X*(n_X+1)/2,H);

A4r = zeros(1,H);
B4r = zeros(n_X,H);
C4r = zeros(n_Z,H);
D4r = zeros(n_X*(n_X+1)/2,H);

a   = zeros(1,1);
b   = zeros(n_X,1);
c   = zeros(n_Z,1);
d   = zeros(n_X*(n_X+1)/2,1);

for h = 1:H
    u_h_1 = [(b - mu_piX);(c - mu_piZ);(d - mu_piXX)];
    PsiQ = compute_LTQ_Y(model_sol,u_h_1);
    %PsiQ = compute_LT_Y(model_sol,u_h_1);
    
    a = - eta0star - mu_pi0 + a + PsiQ.Psi0;
    b = - eta1star + PsiQ.Psi_X;
    c = PsiQ.Psi_Z;
    d = PsiQ.Psi_XX;
    
    A(:,h) = a;
    B(:,h) = b;
    C(:,h) = c;
    D(:,h) = d;
    
    A4r(:,h) = -a/h;
    B4r(:,h) = -b/h;
    C4r(:,h) = -c/h;
    D4r(:,h) = -d/h;
end





