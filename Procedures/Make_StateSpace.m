function[StateSpace] = Make_StateSpace(model_sol,stdv_measur_eq)
% This function casts the model into a state-space representation whose
% format is that that can be used by the "kalman" function.

n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X^2;

nb_measurement_eq = max(size(stdv_measur_eq));

% Compute VAR representation:
H = 1;
[E,V,A,B,Theta0,Theta1,AH,BH,Theta0H,Theta1H] = compute_EV(model_sol,H);

StateSpace.A = zeros(1,nb_measurement_eq);
StateSpace.H = zeros(n_Y,nb_measurement_eq);

StateSpace.mu_xi  = AH;
StateSpace.F      = BH;
StateSpace.Theta0 = Theta0;
StateSpace.Theta1 = Theta1;
StateSpace.R      = diag(stdv_measur_eq.^2);

% Define measurement equations:

nb_measur = 0; % count measurement equations

% Consumption/GDP growth:
nb_measur = nb_measur + 1;
StateSpace.A(1,nb_measur) = model_sol.mu_c0;
vec_dc        = zeros(n_Y,1);
vec_dc(1:n_X) = model_sol.mu_c1;
StateSpace.H(:,nb_measur) = vec_dc;

%Inflation:
nb_measur = nb_measur + 1;
StateSpace.A(1,nb_measur) = model_sol.mu_pi0;
vec_pi  = [model_sol.mu_piX;model_sol.mu_piZ;model_sol.mu_piXX];
StateSpace.H(:,nb_measur) = vec_pi;





