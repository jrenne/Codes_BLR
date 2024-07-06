function[APr,BPr,APrn,BPrn] = compute_AB_real_nom_P(model_sol,H,...
    eta0r,eta1r,eta0rn,eta1rn,...
    A1,B1)
% This functions computes the term structure of rea/nominal rates under P.
% (under the Expectation Hypothesis, without convexity adjustment).

% Real rate:
% eta0star = model_sol.eta0star;
% eta1star = model_sol.eta1star;

n_X = size(model_sol.PhiQ,1);
n_Z = size(model_sol.Phi_Z,1);
n_Y = n_X + n_Z + n_X * (n_X + 1)/2;

I_PhiY_1     = (eye(n_Y) - B1)^(-1);

% Results will be stored in:
APr  = zeros(1,H);
BPr  = zeros(n_Y,H);
APrn = zeros(1,H);
BPrn = zeros(n_Y,H);

APr(1)    = eta0r;
BPr(:,1)  = eta1r;
APrn(1)   = eta0rn;
BPrn(:,1) = eta1rn;

% Compute conditional moments, means across horizons 1 to H ---------------
A = zeros(n_Y,1);
B = eye(n_Y);
for h = 1:(H-1)
    Ah = I_PhiY_1 * (eye(n_Y) - B1^h) * A1;
    Bh = B1^h;
    A  = A + Ah;
    B  = B + Bh;

    %Save results:
    APr(h+1)   = eta0r + (1/(h+1)) * (A' * eta1r);
    BPr(:,h+1) = (1/(h+1)) * B' * eta1r;

    APrn(h+1)   = eta0rn + (1/(h+1)) * (A' * eta1rn);
    BPrn(:,h+1) = (1/(h+1)) * B' * eta1rn;
end

