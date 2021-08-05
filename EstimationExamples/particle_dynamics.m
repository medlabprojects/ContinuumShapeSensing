clear;

% We consider a particle being driven by a sinusoidal forcing function 
m = 1;
f = @(t) sin(t);
T = 10;

% The state is defined by x = [p; p_dot]
state_deriv = @(t, x) [0, 1; 0, 0]*x + [1; 1/m]*f(t);

x_0 = [0; 0];

t = [0:0.1:T];

[t, x_truth] = ode45(state_deriv, t, x_0);

ds = 10;
sigma_y = 1;
t_meas = downsample(t, ds);
y = downsample(x_truth(:,1), ds) + 0.1*randn(size(t_meas));

figure(1);
clf;
hold on;

plot(t, x_truth(:,1));
plot(t_meas, y, '.k', 'MarkerSize', 10);

legend('Truth', 'Data');

x_weg = reshape(downsample(x_truth, ds)', 1, []);