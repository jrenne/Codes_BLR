function[xi_constrained] = impose_constraint(xi,option)

n_X = option;
n_Y = size(xi,1) - n_X - n_X*(n_X+1)/2;

xi_constrained = xi;
x = xi(1:n_X);
xi_constrained((n_X+n_Y+1):end) = vech(x * x');
