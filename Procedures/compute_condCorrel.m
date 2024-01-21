function [coVar,condCorrel] = compute_condCorrel(model_sol,H,vec1,vec2,Y)

T = size(Y,1);
n_Y = size(Y,2);

[~,~,~,~,~,~,~,~,Theta0H,Theta1H] = compute_EV(model_sol,H);

% Compute conditional covariance matrices:
condCov = ones(T,1) * Theta0H' + Y * Theta1H';

coVar = condCov .* (ones(T,1) * kron(vec1',vec2'));
coVar = coVar * ones(n_Y^2,1);

var1 = condCov .* (ones(T,1) * kron(vec1',vec1'));
var1 = var1 * ones(n_Y^2,1);

var2 = condCov .* (ones(T,1) * kron(vec2',vec2'));
var2 = var2 * ones(n_Y^2,1);

condCorrel = coVar ./ (sqrt(var1) .* sqrt(var2));


