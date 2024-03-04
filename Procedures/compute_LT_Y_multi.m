function[Psi] = compute_LT_Y_multi(model,u,H)
% Computes the physical Laplace transform of Y
Psi0 = 0;
v = u;
for i = 1:H
    [Psi] = compute_LT_Y(model,v);
    v     = Psi.Psi_Y;
    Psi0  = Psi0 + Psi.Psi0;
end
Psi.Psi0    = Psi0;
