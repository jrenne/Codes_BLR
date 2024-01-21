function[StateSpace,xi_00,P_00,A4r_raw,B4r_raw,...
    A4rn_raw,B4rn_raw,C4rn_raw,D4rn_raw,...
    eta0r,eta1r,eta0rn,eta1rn] = prepare_State_Space(model_sol,Data_StateSpace)
% This function prepares the state space based on a solved model
% (model_sol) and the dataset

global frequency;

n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

% Compute (conditional) moments of Y:
H = 1; % horizon: 10 years (for SPF forecasts)
[xi_00,P_00,A1,B1,Theta0,Theta1,AH,BH,Theta0H,Theta1H] = compute_EV(model_sol,H);

StateSpace = struct;

StateSpace.mu_xi  = A1;
StateSpace.F      = B1;
StateSpace.Theta0 = Theta0H;
StateSpace.Theta1 = Theta1H;

% Compute interest-rate loadings:
maturities_in_year = Data_StateSpace.maturties_nomyields_in_years;
maturities = maturities_in_year * frequency; % for nominal yields
maturities_real_in_year = Data_StateSpace.maturties_reayields_in_years;
maturities_real = maturities_real_in_year * frequency; % for real yields
H = max([maturities;maturities_real]);
[~,~,A4r,B4r] = compute_AB(model_sol,H);
[~,~,~,~,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,H);

loadings4r  = frequency * [B4r;0*C4rn;0*D4rn];
loadings4rn = frequency * [B4rn;C4rn;D4rn];

A4r_raw  = A4r;
B4r_raw  = B4r;
A4rn_raw = A4rn;
B4rn_raw = B4rn;
C4rn_raw = C4rn;
D4rn_raw = D4rn;

A4r    = frequency * A4r;
A4rn   = frequency * A4rn;

% Save specification of short-term real and nominal rates:
eta0r  = A4r(1)/frequency;
eta1r  = loadings4r(:,1)/frequency;
eta0rn = A4rn(1)/frequency;
eta1rn = loadings4rn(:,1)/frequency;

%     'GDP growth'
vec_dc        = zeros(n_Y,1);
vec_dc(1:n_X) = model_sol.mu_c1;

%     'output gap'
vec_z    = zeros(n_Y,1);
vec_z(2) = 1;

%     'inflation'
vec_pi  = [model_sol.mu_piX;model_sol.mu_piZ;model_sol.mu_piXX];

% ---- SPFs ----
I_PhiY_1     = (eye(n_Y) - B1)^(-1);
H = frequency * 10; % horizon of SPFs
% Compute conditional moments, means across horizons 1 to H ---------------
Amean = zeros(n_Y,1);
Bmean = zeros(n_Y,n_Y);
for h = 1:H
    Ah = I_PhiY_1 * (eye(n_Y) - B1^h) * A1;
    Bh = B1^h;
    Amean = Amean + 1/H * Ah;
    Bmean = Bmean + 1/H * Bh;
end
%  E(inflation 10 years) [measurement variable: Philly Fed SPF]
vec_cpi_10y  = Bmean' * vec_pi ;
cons_cpi_10y = model_sol.mu_pi0 + Amean' * vec_pi ;
%  E(GDP growth 10 years) [measurement variable: Philly Fed SPF]
vec_GDP_10y  = Bmean' * vec_dc ;
cons_GDP_10y = model_sol.mu_c0 + Amean' * vec_dc ;
%  E(BOND 10 years) [measurement variable: Philly Fed SPF]
maturity_in_years_BOND10 = 10 ;
vec_BOND10  = Bmean' * loadings4rn(:,frequency*maturity_in_years_BOND10) ;
cons_BOND10 = A4rn(frequency*maturity_in_years_BOND10) + Amean' * vec_BOND10 ;
%  E(BILL10 years) [measurement variable: Philly Fed SPF]
maturity_in_years_BILL10 = .25;
vec_BILL10  = Bmean' * loadings4rn(:,frequency*maturity_in_years_BILL10) ;
cons_BILL10 = A4rn(frequency*maturity_in_years_BILL10) + Amean' * vec_BILL10 ;


StateSpace.R      = diag(Data_StateSpace.stdv_measur.^2);

StateSpace.A      = [100*model_sol.mu_c0 100*model_sol.mu_pi0 0 ...
    100*A4rn(maturities) ...
    100*A4r(maturities_real) ...
    100 * frequency * cons_cpi_10y ...%100 * frequency * cons_GDP_10y ...
    100 * cons_BOND10,...
    100 * cons_BILL10];
StateSpace.H      = [100*vec_dc 100*vec_pi vec_z ...
    100*loadings4rn(:,maturities) ...
    100*loadings4r(:,maturities_real) ...
    100 * frequency * vec_cpi_10y ...%100 * frequency * vec_GDP_10y ...
    100 * vec_BOND10,...
    100 * vec_BILL10];

StateSpace.option = n_X;






