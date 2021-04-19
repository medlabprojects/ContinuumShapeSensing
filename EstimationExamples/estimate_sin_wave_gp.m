% In this code, we do an example where we estimate a function using GP
% regression. We take samples from a sin wave as ground truth and use
% GP regresssion to estimate the mean and covariance of the function at
% many test points.

clear;

sigma_eps = 0.1;
alpha = 1;

N = 100;
a = 1;
b = 9;

x = linspace(a, b, N)';
x_star = linspace(a - 1, b + 1, 100*N)';

f = sin(x) + sigma_eps*randn(size(x));

f_star_bar = exp_cov(x_star, x, alpha)*inv(exp_cov(x, x, alpha) + sigma_eps^2*eye(length(x)))*f;
f_star_cov = exp_cov(x_star, x_star, alpha) - exp_cov(x_star, x, alpha)*inv(exp_cov(x, x, alpha) + sigma_eps^2*eye(length(x)))*exp_cov(x, x_star, alpha);

figure(1);
clf;
hold on;

plot(x_star, sin(x_star));
plot(x, f, '.k', 'MarkerSize', 10);
plot(x_star, f_star_bar);

p = patch([x_star; flipud(x_star)], [f_star_bar + 2.*sqrt(abs(diag(f_star_cov))); flipud(f_star_bar - 2.*sqrt(abs(diag(f_star_cov))))], [1, 0, 0]);
p.FaceAlpha = 0.1;

legend('Truth', 'Data', 'Mean', 'Covariance');


function K = exp_cov(x_1, x_2, alpha)
    K = zeros(length(x_1), length(x_2));
    
    for iii = 1:length(x_1)
        for jjj = 1:length(x_2)
            A = -1/2*(abs(x_1(iii) - x_2(jjj))^2);
            K(iii,jjj) = alpha*exp(A);
        end
    end
end