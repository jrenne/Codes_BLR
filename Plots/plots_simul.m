
n_X   = size(model_sol.Phi,1);
n_Z   = size(model_sol.Phi_Z,1);
n_eps = size(model_sol.Sigma,2);

% Make sure model solution is updated:
model_sol = compute_sdf(model_sol);

% Simulate model:
T = 200;
x0 = zeros(n_X,1);
z0 = zeros(n_Z,1);
[Delta_c,gamma,inflation,X,Z,XX,Y] = simul_model(model_sol,T,x0,z0);

% Compute (un)conditional moments:
H = 1;
[E,V,A,B,Theta0,Theta1,AH,BH,Theta0H,Theta1H] = compute_EV(model_sol,H);


% Compare simulated and sample variances:
% disp([var(Y)' diag(V)]);

count = 0;

% Plot consumption growth and inflation:
count = count + 1;
subplot(3,3,count);
plot([Delta_c inflation]);
title('Consumption growth and inflation (red)')

% Plot risk aversion:
count = count + 1;
subplot(3,3,count);
plot([gamma,X(:,6)]);
title('Risk Aversion (factor w in red)')

% Plot output gap and potential growth:
count = count + 1;
subplot(3,3,count);
plot([X(:,1) X(:,2) X(:,3)]);
title('Potential growth g and output gap components (d in red, s in orange)')

% Create model where P = Q:
model_P = model_sol;
model_P.muQ = 0 * model_P.muQ;
model_P.PhiQ = model_P.Phi;
model_P.muZQ = 0*model_sol.muZQ;
model_P.Phi_ZXQ = 0*model_sol.Phi_ZXQ;
model_P.Phi_ZXXQ = 0*model_sol.Phi_ZXXQ;

% Compute real and nominal yield curves:
H = 20 * frequency;
[A,B,A4r,B4r] = compute_AB(model_sol,H);
[A_P,B_P,A4r_P,B4r_P] = compute_AB(model_P,H);
[An,Bn,Cn,Dn,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,H);
[An_P,Bn_P,Cn_P,Dn_P,A4rn_P,B4rn_P,C4rn_P,D4rn_P] =...
    compute_AB_nominal(model_P,H);

% Compute simulated real yields:
simr = frequency * (ones(T,1)*A4r + X * B4r);
simr_P = frequency * (ones(T,1)*A4r_P + X * B4r_P);

% Compute simulated nominal yields:
simrn = frequency * (ones(T,1)*A4rn + X * B4rn + Z * C4rn + XX * D4rn);
simrn_P = frequency * (ones(T,1)*A4rn_P + X * B4rn_P + Z * C4rn_P + XX * D4rn_P);

% Plot average TS of real and nominal yields:
count = count + 1;
subplot(3,3,count);
plot((1:H)/frequency,frequency * [(A4r'+B4r'*E(1:n_X)) (A4rn'+[B4rn' C4rn' D4rn']*E)]);
title('Average term structures (real in blue)')

% Plot time series of real yields:
count = count + 1;
subplot(3,3,count);
plot(simr(:,[3,10]));
title('Real yields (long-term in red)')

% Plot time series of real yields (under P):
count = count + 1;
subplot(3,3,count);
plot(simr_P(:,[3,10]));
title('Real yields under P (long-term in red)')

% Plot real term premiums:
count = count + 1;
subplot(3,3,count);
plot(simr(:,[3,10]) - simr_P(:,[3,10]));
title('Real term premiums (long-term in red)')

% Plot nominal term premiums:
count = count + 1;
subplot(3,3,count);
plot(simrn(:,[3,10]) - simrn_P(:,[3,10]));
title('Nominal term premiums (long-term in red)')

% Plot k factor:
count = count + 1;
subplot(3,3,count);
plot(X(:,7));
title('k factor')


