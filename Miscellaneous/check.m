
t = 243;
x0 = all_xi_tt(t,1:n_X)';
z0 = all_xi_tt(t,(n_X+1):(n_X+n_Z))';

plot(all_xi_tt(:,5));

T = 20;
nb_replic = 100000;

[X,Z,XX,Y] = simul_model_multi(model_sol,T+1,x0,z0,nb_replic);

vec_z    = zeros(n_Y,1);
vec_z(2) = 1;
vec_pi_tilde  = zeros(n_Y,1);
vec_pi_tilde(n_X+2) = 1;

z = X(2,:,T+1)';
pi_tilde = Z(2,:,T+1)';

%plot([z pi_tilde]);

disp(corrcoef([z pi_tilde]));
disp(100*cov([z pi_tilde]));

[Theta0H,Theta1H] = compute_condCov(model_sol,T);
% Compute conditional covariance matrices:

Y = all_xi_tt(t,:)';

condCov = reshape(Theta0H + Theta1H * Y,n_Y,n_Y);


A = [vec_z vec_pi_tilde];
disp(100*(A' * condCov * A));
