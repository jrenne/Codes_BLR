function[uncond_r,uncond_rn,uncond_TPr,uncond_TPrn,uncond_IRP,...
    all_stdv_r,all_stdv_rn,all_stdv_TPr,all_stdv_TPrn,all_stdv_IRP] = compute_uncond_yds_TP(model_sol,H)
% Computes conditional and unconditional moments (expectations and variances)

global frequency;

% Factor moments:
[E,V,A1,B1] = compute_EV(model_sol,1); % needed to compute expectations

% Compute uncond yield curves:
[~,~,A4r,B4r] = compute_AB(model_sol,H);
[~,~,~,~,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,H);

loadings4r  = 100 * frequency * [B4r;0*C4rn;0*D4rn];
loadings4rn = 100 * frequency * [B4rn;C4rn;D4rn];
constant4r  = 100 * frequency * A4r;
constant4rn = 100 * frequency * A4rn;

% Save specification of short-term real and nominal rates:
eta0r  = A4r(1);
eta0rn = A4rn(1);
eta1r  = loadings4r(:,1)/(100 * frequency);
eta1rn = loadings4rn(:,1)/(100 * frequency);

% Under Expectation Hypothesis:
[APr,BPr,APrn,BPrn] = compute_AB_real_nom_P(model_sol,H,...
    eta0r,eta1r,eta0rn,eta1rn,A1,B1);
loadings4r_P  = 100 * frequency * BPr;
loadings4rn_P = 100 * frequency * BPrn;
constant4r_P  = 100 * frequency * APr;
constant4rn_P = 100 * frequency * APrn;


% Loadings for TP and IRP:
loadingsTPr  = loadings4r - loadings4r_P;
loadingsTPrn = loadings4rn - loadings4rn_P;
loadingsIRP  = loadingsTPrn - loadingsTPr;

constantTPr  = constant4r - constant4r_P;
constantTPrn = constant4rn - constant4rn_P;
constantIRP  = constantTPrn - constantTPr;

% compute stdv of yields:
all_stdv_r  = zeros(H,1);
all_stdv_rn = zeros(H,1);
all_stdv_TPr  = zeros(H,1);
all_stdv_TPrn = zeros(H,1);
all_stdv_IRP  = zeros(H,1);

for h = 1:H
    all_stdv_r(h)  = sqrt(loadings4r(:,h)'  * V * loadings4r(:,h));
    all_stdv_rn(h) = sqrt(loadings4rn(:,h)' * V * loadings4rn(:,h));

    all_stdv_TPr(h)  = sqrt(loadingsTPr(:,h)'  * V * loadingsTPr(:,h));
    all_stdv_TPrn(h) = sqrt(loadingsTPrn(:,h)' * V * loadingsTPrn(:,h));
    all_stdv_IRP(h)  = sqrt(loadingsIRP(:,h)'  * V * loadingsIRP(:,h));
end

uncond_r  = constant4r'  + loadings4r' * E;
uncond_rn = constant4rn' + loadings4rn' * E;
uncond_TPr  = constantTPr'  + loadingsTPr' * E;
uncond_TPrn = constantTPrn' + loadingsTPrn' * E;
uncond_IRP  = constantIRP'  + loadingsIRP' * E;
