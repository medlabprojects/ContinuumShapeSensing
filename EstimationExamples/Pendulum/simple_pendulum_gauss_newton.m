% Script for simulation and then performing continuous-time batch
% estimation for the case of a simple pendulum whose angle coordinate is
% observed with an unknown offset. Given pendulum inputs, 
% timestamped measurements for the horizontal distance, we want to estimate
% the state of the pendulum vs. time as well as the angulur offset of the 
% endocoder.

clear;

% First we simulate the pendulum motion
u = @(t) 0;
Q_mag = 1e1;
Q = diag([Q_mag^2, Q_mag^2]);
x_0 = [2; 0];
t_0 = 0;
t_f = 15;
t_span = [t_0, t_f];
T = 0.001;
m_gt = 0.25;
[t_gt, x_gt, z_gt] = simulate_pendulum_motion(x_0, u, t_span, T, m_gt, Q);

% Now subsample the data
ds = 100;
z = downsample(z_gt, ds);
t = downsample(t_gt, ds);

R_z = (5e-2)^2;

% Add noise to the data
z = z + mvnrnd(0, R_z, length(z))';

% Now do the estimation using gauss newton
spline = MatrixSpline();
spline.setup(t_0, t_f, 50, 2);

max_iterations = 30;
del_theta_tol = 1e-8;

m_hat = 0;
P_m = (0.5)^2;
x_0_hat = [2.75; 0];
P_x = diag([5^2, 5^2]);
Q = diag([2e-1^2, 2e-1^2]);

theta_bar = [spline.coefficients; m_hat];

for iii = 1:max_iterations
    c_bar = theta_bar(1:(end - 1));
    spline.set_coefficients(c_bar);
    
    [A, b] = build_A_b(theta_bar, spline, t, z, R_z, m_hat, P_m, t_0, t_f, x_0_hat, P_x, Q, u);
    del_theta = A \ (-b);
    theta_bar = theta_bar + del_theta;
    
    fprintf('norm(del_theta): %f\n', norm(del_theta));
    
    if norm(del_theta) < del_theta_tol
        break;
    end
end

theta_star = theta_bar;
theta_cov = inv(A);
c_cov = theta_cov(1:(end - 1), 1:(end - 1));

x_cov = @(t) spline.get_phi(t, 0)*c_cov*(spline.get_phi(t, 0)');

x_std = zeros(2,length(t_gt));

for iii = 1:length(t_gt)
    x_std_i = sqrt(diag(x_cov(t_gt(iii))));
    x_std(:,iii) = x_std_i;
end

c_star = theta_star(1:(end - 1));
m_star = theta_star(end);

spline.set_coefficients(c_star);
x_star = spline.evaluate(t_gt, 0);

x_conf_hi = x_star + 2.*x_std;
x_conf_lo = x_star - 2.*x_std;

theta_conf_lo = x_conf_lo(1,:);
theta_conf_hi = x_conf_hi(1,:);
theta_dot_conf_lo = x_conf_lo(2,:);
theta_dot_conf_hi = x_conf_hi(2,:);

t_conf = [t_gt, t_gt(end:-1:1)];
theta_conf = [theta_conf_hi, theta_conf_lo(end:-1:1)];
theta_dot_conf = [theta_dot_conf_hi, theta_dot_conf_lo(end:-1:1)];

z_star = zeros(1,length(t_gt));

for iii = 1:length(t_gt)
    z_star(iii) = observation_equation(x_star(:,iii), m_star);
end

figure(1);
clf;
hold on;

plot(t_gt, x_gt);
plot(t_gt, x_star);

fill(t_conf, theta_conf, 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill(t_conf, theta_dot_conf, 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

legend('theta truth', 'thetaDot truth', 'theta star', 'thetaDot star');

figure(2);
clf;
hold on;

plot(t, z, '.k', 'MarkerSize', 10);
plot(t_gt, z_gt);
plot(t_gt, z_star);

legend('z', 'zGt', 'zStar');


function [A, b] = build_A_b(theta, spline, t, z, R_z, m_hat, P_m, t_0, t_f, x_0_hat, P_x, Q, u)
    [A_z, b_z] = build_A_b_z(theta, t, z, R_z, spline);
    [A_m, b_m] = build_A_b_m(theta, m_hat, P_m);
    [A_x, b_x] = build_A_b_x(theta, spline, t_0, x_0_hat, P_x);
    [A_u, b_u] = build_A_b_u(theta, spline, Q, u, t_0, t_f);
    
    A = sum(A_z, 3) + A_m + A_x + A_u;
    b = sum(b_z, 2) + b_m + b_x + b_u;
end

function [A_z, b_z] = build_A_b_z(theta, t, z, R_z, spline)
    % Now, there'll be matrices for each measurement, so these will be 3D
    A_z = zeros(length(theta), length(theta), size(z,2));
    b_z = zeros(length(theta), size(z,2));
    
    m = theta(end);
    
    for iii = 1:length(z)
        [J_x, J_m] = observation_jacobian(spline.evaluate(t(iii), 0), m);
        
        phi_i = spline.get_phi(t(iii),0);
        E_z_i = [-J_x*phi_i, -J_m];
        e_z_bar_i = z(:,iii) - observation_equation(spline.evaluate(t(iii),0), m);
        
        A_z(:,:,iii) = (E_z_i')*inv(R_z)*E_z_i;
        b_z(:,iii) = (E_z_i')*inv(R_z)*e_z_bar_i;
    end
end

function [A_m, b_m] = build_A_b_m(theta, m_hat, P_m)
    c = theta(1:(end - 1));
    m = theta(end);
    
    E_m = [zeros(length(m), length(c)), eye(length(m))];
    e_m_bar = m - m_hat;
    
    A_m = (E_m')*inv(P_m)*E_m;
    b_m = (E_m')*inv(P_m)*e_m_bar;
end

function [A_x, b_x] = build_A_b_x(theta, spline, t_0, x_0_hat, P_x)
    m = theta(end);
    
    phi_0 = spline.get_phi(t_0, 0);
    
    E_x = [phi_0, zeros(size(phi_0,1),length(m))];
    e_x_bar = spline.evaluate(t_0, 0) - x_0_hat;
    
    A_x = (E_x')*inv(P_x)*E_x;
    b_x = (E_x')*inv(P_x)*e_x_bar;
end

function [A_u, b_u] = build_A_b_u(theta, spline, Q, u, t_0, t_f)
    
    m = theta(end);
    
    J = @(t) state_jacobian(spline.evaluate(t,0), u(t));
    
    E_u = @(t) [spline.get_phi(t, 1) - J(t)*spline.get_phi(t, 0), zeros(2, length(m))];
    e_u_bar = @(t) spline.evaluate(t, 1) - state_equation(spline.evaluate(t, 0), u(t));
    
    % Now we need to integrate two things to get A and B
    integrand_A = @(t) (E_u(t)')*inv(Q)*E_u(t);
    integrand_b = @(t) (E_u(t)')*inv(Q)*e_u_bar(t);
    
    ttt = t_0:0.01:t_f;
    AAA = zeros(length(theta), length(theta), length(ttt));
    bbb = zeros(length(theta), length(ttt));
    
    for iii = 1:length(ttt)
        AAA(:,:,iii) = integrand_A(ttt(iii));
        bbb(:,iii) = integrand_b(ttt(iii));
    end
    
    A_u = trapz(ttt, AAA, 3);
    b_u = trapz(ttt, bbb, 2);
end

function J = finite_difference_jacobian(f, x)
    del = 1e-6;
    
    f_0 = f(x);
    J = zeros(length(f_0), length(x));
    
    for iii = 1:length(x)
        del_x = zeros(size(x));
        del_x(iii) = del;
        del_f = f(x + del_x) - f_0;
        J(:,iii) = del_f/del;
    end
end

function [J_x, J_m] = observation_jacobian(x, m)
    f = @(x) observation_equation(x, m);
    J_x = finite_difference_jacobian(f, x);
    
    f = @(m) observation_equation(x, m);
    J_m = finite_difference_jacobian(f, m);
end

function J = state_jacobian(x, u)
    f = @(x) state_equation(x, u);
    J = finite_difference_jacobian(f, x);
end