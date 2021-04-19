% This code is an implementation of the GP sampling outline given by eq.
% 2-17 of Rassmussen and Williams.

% First pick where we want to sample.
x_star = 0:0.01:10;

% Next build the kernel (covariance matrix of the samples)
K = exp_covariance(x_star, x_star);

% Sample from the gaussian specified by 2-17
num_samples = 5;
f_star = mvnrnd(zeros(length(x_star),1)', K, num_samples);

% Plot the results
figure(1);
clf;
hold on;

plot(x_star, f_star);



function K = exp_covariance(x_1, x_2)
    K = zeros(length(x_1), length(x_2));
    
    for iii = 1:length(x_1)
        for jjj = 1:length(x_2)
            A = -1/2*(abs(x_1(iii) - x_2(jjj))^2);
            K(iii,jjj) = exp(A);
        end
    end
end