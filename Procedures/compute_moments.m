function[mim] = compute_moments(model_sol,A4r,B4r,A4rn,B4rn,C4rn,D4rn)
% This functions computes the distance between model-implied moments and targets.

global frequency;

n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

% Compute moments of Y:
H = 1;
[E,V,A1,B1,~,~,~,BH,~,~] = compute_EV(model_sol,H);

% Compute covariance matrix of [Y_t;Y_{t+1}]:
VarYttp1 = V * BH';

% Compute interest rate loadings:
maturities_in_year = [2;10];
maturities = maturities_in_year * frequency;
H = max(maturities);
%[~,~,A4r,B4r] = compute_AB(model_sol,H);
%[~,~,~,~,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,H);
loadings4r  = frequency * [B4r;0*C4rn;0*D4rn];
loadings4rn = frequency * [B4rn;C4rn;D4rn];

A4r    = frequency * A4r;
A4rn   = frequency * A4rn;

% Average real rates:
avg_r       = A4r + E' * loadings4r;

% Average nominal rates:
avg_rn      = A4rn + E' * loadings4rn;

% Compute loadings under EH:
loadings4r  = frequency * [B4r;0*C4rn;0*D4rn];
loadings4rn = frequency * [B4rn;C4rn;D4rn];

% Keep specification of short-term real and nominal rates:
eta0r  = A4r(1)/frequency;
eta1r  = loadings4r(:,1)/frequency;
eta0rn = A4rn(1)/frequency;
eta1rn = loadings4rn(:,1)/frequency;
[~,~,~,BPrn] = compute_AB_real_nom_P(model_sol,H,...
    eta0r,eta1r,eta0rn,eta1rn,A1,B1);
% Compute loadings for nominal rates under EH:
loadings4rn_P = frequency * BPrn;

% Variance of real yields:
var_r = loadings4r' * V * loadings4r;

% Variance of nominal yields:
var_rn = loadings4rn' * V * loadings4rn;

% % Compute correlation between risk aversion and consumption growth:
% cov_RA_deltaC = model_sol.mu_c1' * V(1:n_X,1:n_X) * model_sol.mu_gamma1;
% cor_RA_deltaC = cov_RA_deltaC / (stdv_dc * stdv_gamma);

% Construct vector of model-implied moments:
mim = [];

% Avg GDP growth
mim           = [mim;model_sol.mu_c0];

% Avg inflation
mim           = [mim;model_sol.mu_pi0];

%     'Std. dev. of GDP growth';...
vec_dc        = zeros(n_Y,1);
vec_dc(1:n_X) = model_sol.mu_c1;
stdv_dc       = sqrt(vec_dc' * V * vec_dc);
mim           = [mim;stdv_dc];

%     'Std. dev. of output gap';...
vec_z    = zeros(n_Y,1);
vec_z(2) = 1;
stdv_z   = sqrt(vec_z' * V * vec_z);
mim      = [mim;stdv_z];

%     'Std. dev. of inflation';...
vec_pi  = [model_sol.mu_piX;model_sol.mu_piZ;model_sol.mu_piXX];
stdv_pi = sqrt(vec_pi' * V * vec_pi);
mim     = [mim;stdv_pi];

%     'AC(1) of GDP growth';...
covar  = vec_dc' * VarYttp1 * vec_dc;
AC1_dc = covar/stdv_dc^2;
mim = [mim;AC1_dc];

%     'AC(1) of output gap';...
covar   = vec_z' * VarYttp1 * vec_z;
AC1_z   = covar/stdv_z^2;
mim = [mim;AC1_z];

%     'AC(1) of inflation';...
covar   = vec_pi' * VarYttp1 * vec_pi;
AC1_pi  = covar/stdv_pi^2;
mim = [mim;AC1_pi];

%     'Min of ratio E(RA)/SD(RA)';
vec_gamma        = zeros(n_Y,1);
vec_gamma(1:n_X) = model_sol.mu_gamma1;
stdv_gamma       = sqrt(vec_gamma' * V * vec_gamma);
mim = [mim;model_sol.mu_gamma0/stdv_gamma];

%     'Max of E(RA)';...
mim = [mim;model_sol.mu_gamma0];

%     'Max correlation between output gap and RA';...
covar   = vec_z' * VarYttp1 * vec_gamma;
mim = [mim;covar/(stdv_z*stdv_gamma)];

%     'Avg 10-year real yield';...
mim = [mim;avg_r(maturities(2))];

%     'Avg 3-month nominal yield';...
mim = [mim;avg_rn(1)];

%     'Avg 10-year nominal yield';...
mim = [mim;avg_rn(maturities(2))];

%     'Std. dev. of 10-year real rate';...
var_ltr = loadings4r(:,maturities(2))' * V * loadings4r(:,maturities(2));
mim = [mim;sqrt(var_ltr)];

%     'Std. dev. of 10-year nominal rate';...
var_ltrn = loadings4rn(:,maturities(2))' * V * loadings4rn(:,maturities(2));
mim  = [mim;sqrt(var_ltrn)];

%     'Avg slope of real yield curve (2-10 yrs)';
mim = [mim;avg_r(maturities(2))-avg_r(maturities(1))];

%     'Avg slope of nominal yield curve (2-10 yrs)';...
mim = [mim;avg_rn(maturities(2))-avg_rn(maturities(1))];

%     'Std. dev. of 2-year real rate';...
var_str = loadings4r(:,maturities(1))' * V * loadings4r(:,maturities(1));
mim = [mim;sqrt(var_str)];

%     'Std. dev. of 2-year nominal rate';...
var_strn = loadings4rn(:,maturities(1))' * V * loadings4rn(:,maturities(1));
mim  = [mim;sqrt(var_strn)];

%     'AC(1) of 10-year real yield';...
covar   = loadings4r(:,maturities(2))' * VarYttp1 * loadings4r(:,maturities(2));
mim = [mim;covar/var_ltr];

%     'AC(1) of 10-year nominal rate';...
covar    = loadings4rn(:,maturities(2))' * VarYttp1 * loadings4rn(:,maturities(2));
mim  = [mim;covar/var_ltrn];

%     'Std. dev. of 10-year nominal term premium';...
vec_tpn  = loadings4rn(:,maturities(2))-loadings4rn_P(:,maturities(2));
stdv_tpn = sqrt(vec_tpn' * V * vec_tpn);
mim = [mim;stdv_tpn];

%     'Correl. 10-year real rate and output gap'
covar   = loadings4r(:,maturities(2))' * V * vec_z;
mim = [mim;covar/(stdv_z*sqrt(var_ltr))];

