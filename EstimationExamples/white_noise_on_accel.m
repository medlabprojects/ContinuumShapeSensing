clear;

% The goal here is to simulate the LTI system p_ddot = w(t). What I want to
% do here is:
%   1. Determine the priors on x
%   2. Simulate the system for a specific white noise profile which will
%      give me ground truth data for estimation.
%   3. Solve for the optimal x at the measurement times.

% We use the state variables p and p_dot. With these, the system becomes
% the good old fashioned:
%   x_dot = A*x + B*u + L*w where

A = [0, 1; 0, 0];
B = 0;
L = [0; 1];

% We need this function for sure
exp_A = @(del_t) [1, del_t; 0, 1];
phi = @(t, s) [1, (t - s); 0, 1];

Q = 0.1;
Q_k = @(dt) Q*[1/3*dt^3, 1/2*dt^2; 1/2*dt^2, dt];

% Let's say we measure the thing at these time points
t_K = 1:6;

% Build the large matrix A
A_full = [eye(2),              zeros(2),            zeros(2),            zeros(2),            zeros(2),            zeros(2)
          phi(t_K(2), t_K(1)), eye(2),              zeros(2),            zeros(2),            zeros(2),            zeros(2)
          phi(t_K(3), t_K(1)), phi(t_K(3), t_K(2)), eye(2),              zeros(2),            zeros(2),            zeros(2)
          phi(t_K(4), t_K(1)), phi(t_K(4), t_K(2)), phi(t_K(4), t_K(3)), eye(2),              zeros(2),            zeros(2)
          phi(t_K(5), t_K(1)), phi(t_K(5), t_K(2)), phi(t_K(5), t_K(3)), phi(t_K(5), t_K(4)), eye(2),              zeros(2)
          phi(t_K(6), t_K(1)), phi(t_K(6), t_K(2)), phi(t_K(6), t_K(3)), phi(t_K(6), t_K(4)), phi(t_K(6), t_K(5)), eye(2)];
      
% Let's say that x_0_wed = [0; 0];
v = zeros(size(A_full,1),1);
x_wed = A_full*v;

% Now we also need to build the full Q. Let
P_wed_0 = eye(2);
Q_full = blkdiag(P_wed_0, Q_k(1), Q_k(1), Q_k(1), Q_k(1), Q_k(1));
P_wed = A_full*Q_full*(A_full');

p_wed = x_wed(1:2:end);
tmp = diag(P_wed);
p_wed_var = tmp(1:2:end);

% Our measurement equation is y_k = C_k*x + n_k
R = 0.01;
R_full = R*eye(6);

C = [1, 0];
C_full = blkdiag(C, C, C, C, C, C);

% Now our truth model can actually be anything. Let's choose it as a sin
% wave just for fun
y = sin(t_K)';

P_hat = inv(inv(P_wed) + (C_full')*inv(R_full)*C_full);
x_hat = P_hat*(inv(P_wed)*x_wed + (C_full')*inv(R_full)*y);

p_hat = x_hat(1:2:end);
tmp = diag(P_hat);
p_hat_var = tmp(1:2:end);
p_hat_std = sqrt(p_hat_var);

figure(1);
clf;
hold on;

plot(t_K, p_hat);
xxx = [t_K, fliplr(t_K)];
yyy = [p_hat + p_hat_std; flipud(p_hat - p_hat_std)];
patch(xxx, yyy, 'r', 'FaceAlpha', 0.2);
