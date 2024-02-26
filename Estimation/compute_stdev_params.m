% =========================================================================
% Compute standard deviation of parameter estimates
% =========================================================================

FILTER = FILTER_MLE;

% Create vector of parameters:
parameters = model_sol.param;

Param2remove = (abs(parameters)>10)&(FILTER);
FILTER = FILTER - Param2remove;

sub_parameters =  model_sol.param(FILTER==1);

f = @(x)compute_logl(x,Data_StateSpace,model_sol,indic_add_moments);

% Check log-lik.:
f0 = f(sub_parameters);
disp(f0);

epsilon = .01;

J = zeros(sum(FILTER),sum(FILTER));

for i=1:sum(FILTER)
    disp(i);
    sub_parameters_i = sub_parameters;
    sub_parameters_i(i) = sub_parameters_i(i) + epsilon;
    fi = f(sub_parameters_i);
    for j=1:i
        sub_parameters_ij = sub_parameters_i;
        sub_parameters_ij(j) = sub_parameters_ij(j) + epsilon;
        fij = f(sub_parameters_ij);
        %list_param = find(cumsum(FILTER)==i);
        %indic_param = list_param(1);
        %disp(indic_param);
        J(i,j) = (fij - 2*fj + f0)/epsilon^2;
        J(j,i) = (fij - 2*fj + f0)/epsilon^2;
    end
end

J_1 = J^(-1);
disp([sub_parameters sqrt(diag(J_1))]);


