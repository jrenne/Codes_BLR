

% Compute uncondtional yield curves:
H = frequency * 10;
[uncond_r,uncond_rn,uncond_TPr,uncond_TPrn,uncond_IRP,...
    all_stdv_r,all_stdv_rn,all_stdv_TPr,all_stdv_TPrn,all_stdv_IRP] =...
    compute_uncond_yds_TP(model_sol,H);

% Check what happens when rho_gz = 0 (no hyteresis)
indic_variable = find(strcmp(model_sol.names_param,'rho_gz'));
model_noHysteresis = model_sol;
model_noHysteresis.param(indic_variable) = 0.00001;
model_noHysteresis_sol = make_model_sol(model_noHysteresis);
[uncond_r_noHyteresis] = compute_uncond_yds_TP(model_noHysteresis_sol,H);

n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

% Compute moments of Y:
H = 200;
%H = 80;
[E,V,A,B,Theta0,Theta1,AH,BH,Theta0H,Theta1H] = compute_EV(model_sol,H);
[E_noH,V_noH,A_noH,B_noH,Theta0_noH,Theta1_noH,AH_noH,...
    BH_noH,Theta0H_noH,Theta1H_noH] = compute_EV(model_noHysteresis_sol,H);

% compute conditional variance of Delta c:
vec_dc        = zeros(n_Y,1);
vec_dc(1:n_X) = model_sol.mu_c1;
stdv_dc       = sqrt(vec_dc' * V * vec_dc);

Vh = reshape(Theta0H+Theta1H*E,n_Y,n_Y);
stdv_dc_h       = sqrt(vec_dc' * Vh * vec_dc);

Vh_noH = reshape(Theta0H_noH+Theta1H_noH*E_noH,n_Y,n_Y);
stdv_dc_h_noH       = sqrt(vec_dc' * Vh_noH * vec_dc);

% Conditional variance
Omega = model_sol.Sigma * model_sol.Sigma';
VcumH = zeros(n_X,n_X);
Phik = eye(n_X);
cumPhik = eye(n_X);
for i = 1:H
    VcumH = VcumH + cumPhik * Omega * cumPhik';
    Phik = Phik * model_sol.Phi;
    cumPhik = cumPhik + Phik;
end

Omega_noH = model_noHysteresis_sol.Sigma * model_noHysteresis_sol.Sigma';
VcumH_noH = zeros(n_X,n_X);
Phik = eye(n_X);
cumPhik = eye(n_X);
for i = 1:H
    VcumH_noH = VcumH_noH + cumPhik * Omega_noH * cumPhik';
    Phik = Phik * model_noHysteresis_sol.Phi;
    cumPhik = cumPhik + Phik;
end

disp([sqrt(VcumH(1,1)) sqrt(VcumH_noH(1,1))]);
disp([sqrt(vec_dc(1:n_X)' * VcumH * vec_dc(1:n_X)) sqrt(vec_dc(1:n_X)' * VcumH_noH * vec_dc(1:n_X))]);

Vh_noH = reshape(Theta0H_noH+Theta1H_noH*E_noH,n_Y,n_Y);
stdv_dc_h_noH       = sqrt(vec_dc' * Vh_noH * vec_dc);


FILTER = 0 * model_sol.param;
FILTER(indic_variable) = 1; % stdv_k, mu_c, mu_pi
% Create vector of parameters:
sub_parameters =  model_sol.param(FILTER==1);
% Computation of distance at initial parameters:
Loss_withMoments    = compute_logl(sub_parameters,Data_StateSpace,model_sol,1);
Loss_withoutMoments = compute_logl(sub_parameters,Data_StateSpace,model_sol,0);

sub_parameters =  model_noHysteresis_sol.param(FILTER==1);
Loss_withMoments_noH    = compute_logl(sub_parameters,Data_StateSpace,model_sol,1);
Loss_withoutMoments_noH = compute_logl(sub_parameters,Data_StateSpace,model_sol,0);


