function [process, X, t] = sample_normal_process(x_cov, t_span, dt)
% Assume zero mean process
t = t_span(1):dt:t_span(2);
X = mvnrnd(zeros(length(x_cov), 1), x_cov, length(t))';

% Fit 5th order spline through the data
d = 5;
process = BSpline.interpolate(X, t, d);
end