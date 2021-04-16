clear;

%% Create the normal process function and samples
x_cov = [1, 0.2; 0.2, 1];
t_span = [0, 1];
dt = 0.005;

[process, X, t] = sample_normal_process(x_cov, t_span, dt);

t_s = linspace(t_span(1), t_span(2), 1e4);
X_s = process.evaluate(t_s);

%% Plot results

X_spline = process.evaluate(t);

figure(1);
clf;
hold on;

plot(t_s, X_s);
plot(t, X, '.k');

legend('x_1(t)', 'x_2(t)', 'x_s');