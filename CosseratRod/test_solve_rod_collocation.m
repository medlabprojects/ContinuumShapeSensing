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
N = 50; % Num collocation points, number of conditions is 3*N
M = 12; % Dimension of spline spaces, number of variables is 3*M
s_span = [0, 1];

% Now we need a forcing function f,l : [0,1] -> R3
s = linspace(0, 1, 10000);

% Want guassian loading functions
f1 = normpdf(s, 0.2, 0.01) + normpdf(s, 0.8, 0.01);
f2 = normpdf(s, 0.5, 0.01);

f_data = zeros(3, length(s));
l_data = zeros(3, length(s));

f_data(1,:) = f1;
f_data(2,:) = f2;

d = 5;
f_spline = BSpline.least_squares_fit(f_data, s, d, 500);
l_spline = BSpline.least_squares_fit(l_data, s, d, 500);

f = @(s) f_spline.evaluate(s);
l = @(s) l_spline.evaluate(s);

% Solve and plot the rod
[p, q, n, m] = solve_rod_collocation(f, l, K, N, M, L, p_0, R_0);

figure(1);
clf;

plot_static_rod(p, q, f, l, L);