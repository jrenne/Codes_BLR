function[loadings] = CPregressions(model_sol,fwd_maturities_in_years,...
    xs_maturities_in_years,period_per_year)
% This functions calculates the population parameters of
% the regressions proposed by Cochrane and Piazzesi (2005)
% the one-year maturity is added to the list of forwards

% Compute (un)conditional moments:
H = 1;
[E,V,~,~,~,~,~,BH] = compute_EV(model_sol,H);

% Compute loadings that will be use to compute forward rates:
H = max([fwd_maturities_in_years xs_maturities_in_years]);
[An,Bn,Cn,Dn,A4rn,B4rn,C4rn,D4rn] = compute_AB_nominal(model_sol,H);

A_fwd = An(:,(fwd_maturities_in_years-1)*period_per_year) - ...
    An(:,fwd_maturities_in_years*period_per_year);
B_fwd = Bn(:,(fwd_maturities_in_years-1)*period_per_year) - ...
    Bn(:,fwd_maturities_in_years*period_per_year);
C_fwd = Cn(:,(fwd_maturities_in_years-1)*period_per_year) - ...
    Cn(:,fwd_maturities_in_years*period_per_year);
D_fwd = Dn(:,(fwd_maturities_in_years-1)*period_per_year) - ...
    Dn(:,fwd_maturities_in_years*period_per_year);

% Add one-year maturity:
A = [-An(:,period_per_year) A_fwd];
B = [-Bn(:,period_per_year) B_fwd];
C = [-Cn(:,period_per_year) C_fwd];
D = [-Dn(:,period_per_year) D_fwd];

% We express the returns/yields/forward as combinations of 
% [X_{t},Z_{t},vecXX_{t},X_{t+1},Z_{t+1},vecXX_{t+1}]

loadings_fwd = [A;B;C;D;0*B;0*C;0*D];

A1y = A4rn(:,period_per_year) * ones(1,size(xs_maturities_in_years,2));
B1y = B4rn(:,period_per_year) * ones(1,size(xs_maturities_in_years,2));
C1y = C4rn(:,period_per_year) * ones(1,size(xs_maturities_in_years,2));
D1y = D4rn(:,period_per_year) * ones(1,size(xs_maturities_in_years,2));

loadings_xs = [An(:,(xs_maturities_in_years-1)*period_per_year) - ...
    An(:,xs_maturities_in_years*period_per_year)-A1y;...
    -Bn(:,xs_maturities_in_years*period_per_year)-B1y;...
    -Cn(:,xs_maturities_in_years*period_per_year)-C1y;...
    -Dn(:,xs_maturities_in_years*period_per_year)-D1y;...
    Bn(:,(xs_maturities_in_years-1)*period_per_year);...
    Cn(:,(xs_maturities_in_years-1)*period_per_year);...
    Dn(:,(xs_maturities_in_years-1)*period_per_year)];

n_Y = size(E,1);

VarXXtp1 = kron(eye(2),V);
VarXXtp1(1:n_Y,(n_Y+1):(2*n_Y)) = V * BH';
VarXXtp1((n_Y+1):(2*n_Y),1:n_Y) = BH * V;

EX = [E;E];
EXpX = VarXXtp1 + EX * EX';
Q = [[1;EX],[EX';EXpX]];

XX_1 = (loadings_fwd' * Q * loadings_fwd)^(-1);
XpY = loadings_fwd' * Q * loadings_xs;

loadings = XX_1 * XpY;

plot(loadings);

