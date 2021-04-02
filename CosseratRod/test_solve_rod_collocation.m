clear;

% Rod Properties
D = 1.5/1000;
E = 2.07e11; % n/m^2
G = 7.93e10; % n/m^2
L = 1; % m

I = pi/4*(D/2)^4; % m^4
K = [E*I 0 0; 0 E*I 0; 0 0 2*G*I];

% Integration constants
p_0 = zeros(3,1);
R_0 = eye(3);
N = 100; % Num collocation points
M = 15; % Dimension of spline spaces
s_span = [0, 1];

% Now we need a forcing function f,l : [0,1] -> R3
s = linspace(0, 1, 10000);
f_data = zeros(3, 10000);
l_data = zeros(3, 10000);

% f_data(1:2,55:65) = 5;
f_data(1:2,5000:6000) = -3;

d = 5;
f_spline = BSpline.least_squares_fit(f_data, s, d, 1000);
l_spline = BSpline.least_squares_fit(l_data, s, d, 1000);

f = @(s) f_spline.evaluate(s);
l = @(s) l_spline.evaluate(s);

% Solve and plot the rod
x = solve_rod_collocation(f, l, K, N, M, L, p_0, R_0);

figure(1);
clf;

plot_static_rod(x, f, l, L);