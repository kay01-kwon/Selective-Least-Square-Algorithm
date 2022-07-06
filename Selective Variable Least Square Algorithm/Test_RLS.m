close all
clear all
N = 1000;
x = ones(N,1);
x_est = ones(N,1);

sigma = 0.1;
x_var = 0.010^2;
x_var_vec = [x_var];


for i = 1:N-1
    K = x_var/(x_var+sigma);
    x_est(i+1) = x(i) + K*(x(i)-x_est(i));
    x_var = (1-K)*x_var;
    x_var_vec = [x_var_vec x_var];
end

plot(sqrt(x_var_vec))