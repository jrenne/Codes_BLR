function[Theta0H,Theta1H] = compute_condCov(model,H)
% Computes conditional and unconditional moments (expectations and variances)

epsilon = 10^(-5); % -6

n_X   = size(model.Phi,1);
n_Z   = size(model.Phi_Z,1);

n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

u   = zeros(n_Y,1);

Theta0   = zeros(n_Y*n_Y,1);
Theta1   = zeros(n_Y*n_Y,n_Y);

for i = 1:size(u,1)
    for j = 1:size(u,1)

        u_eps_ipjp  = u;
        u_eps_imjm  = u;
        u_eps_imjp  = u;
        u_eps_ipjm  = u;

        u_eps_ipjp(i)  =  epsilon;
        u_eps_imjm(i)  = -epsilon;
        u_eps_imjp(i)  = -epsilon;
        u_eps_ipjm(i)  =  epsilon;

        u_eps_ipjp(j)  = u_eps_ipjp(j) + epsilon;
        u_eps_imjm(j)  = u_eps_imjm(j) - epsilon;
        u_eps_imjp(j)  = u_eps_imjp(j) + epsilon;
        u_eps_ipjm(j)  = u_eps_ipjm(j) - epsilon;

        LTY_eps_ipjp  = compute_LT_Y_multi(model,u_eps_ipjp,H);
        LTY_eps_imjm  = compute_LT_Y_multi(model,u_eps_imjm,H);
        LTY_eps_imjp  = compute_LT_Y_multi(model,u_eps_imjp,H);
        LTY_eps_ipjm  = compute_LT_Y_multi(model,u_eps_ipjm,H);

        Theta0((i-1)*n_Y+j,1) = (LTY_eps_ipjp.Psi0 + LTY_eps_imjm.Psi0 - ...
            (LTY_eps_imjp.Psi0 + LTY_eps_ipjm.Psi0))/(4*epsilon^2);
        Theta1((i-1)*n_Y+j,:) = (LTY_eps_ipjp.Psi_Y + LTY_eps_imjm.Psi_Y - ...
            (LTY_eps_imjp.Psi_Y + LTY_eps_ipjm.Psi_Y))/(4*epsilon^2);

    end
end

% Compute conditional moments, horizon H ----------------------------------

Theta0H = Theta0;
Theta1H = Theta1;







