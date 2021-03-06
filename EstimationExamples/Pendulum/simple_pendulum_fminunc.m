% Script for simulation and then performing continuous-time batch
% estimation for the case of a simple pendulum whose angle coordinate is
% observed with an unknown offset. Given pendulum inputs, 
% timestamped measurements for the horizontal distance, we want to estimate
% the state of the pendulum vs. time as well as the angulur offset of the 
% endocoder.

clear;

% First we simulate the pendulum motion
t_0 = 0;
t_f = 5;
T = 0.001;
t_gt = t_0:T:t_f;
x_gt = zeros(2, length(t_gt));
z_gt = zeros(1, length(t_gt));
m_gt = 0.25;

x_0 = [2; 0];
x_gt(:,1) = x_0;
z_gt(:,1) = observation_equation(x_0, m_gt);

Q_mag = 5e0;
Q = diag([Q_mag^2, Q_mag^2]);

u = @(t) 0;

for iii = 2:length(t_gt)
    x_dot = state_equation(x_gt(:,iii - 1), u(t_gt(iii))) + mvnrnd([0;0], Q, 1)';
    x_gt(:,iii) = x_gt(:,iii - 1) + x_dot*T;
    z_gt(:,iii) = observation_equation(x_gt(:,iii), m_gt);
end

% Now subsample the data
ds = 1000;
z = downsample(z_gt, ds);
t = downsample(t_gt, ds);

R_z = (5e-2)^2;

% Add noise to the data
z = z + mvnrnd(0, R_z, length(z))';

% Now we're going to do the fit using a nonlinear optimization

spline = MatrixSpline();
spline.setup(t_0, t_f, 10, 2);

c_0 = spline.coefficients;
m_0 = 0;

P_m = (0.5)^2;
x_0 = [2.75; 0];
P_x = diag([5^2, 5^2]);

theta_0 = [c_0; m_0];

obj = @(theta) objective_function(theta, t, z, R_z, spline, m_0, P_m, x_0, P_x, u, diag([0.001, 0.001]));

opts = optimoptions('fminunc', 'disp', 'iter', 'OptimalityTolerance', 1e-8, 'MaxFunctionEvaluations', 1e5, 'MaxIterations', 1e5);

theta_star = fminunc(obj, theta_0, opts);

c_star = theta_star(1:(end - 1));
m_star = theta_star(end);

spline.set_coefficients(c_star);
x_star = spline.evaluate(t_gt, 0);
z_star = zeros(1,length(t_gt));

for iii = 1:length(t_gt)
    z_star(iii) = observation_equation(x_star(:,iii), m_star);
end

figure(1);
clf;
hold on;

plot(t_gt, x_gt);
plot(t_gt, x_star);

legend('theta truth', 'thetaDot truth', 'theta star', 'thetaDot star');

figure(2);
clf;
hold on;

plot(t, z, '.k', 'MarkerSize', 10);
plot(t_gt, z_gt);
plot(t_gt, z_star);

legend('z', 'zGt', 'zStar');


% Now do the estimation using gauss newton
spline.set_coefficients(zeros(size(spline.coefficients)));

max_iterations = 20;

theta_bar = theta_0;

for iii = 1:max_iterations
    [A, b] = build_A_b(theta_bar);
    del_theta = A \ (-b);
    theta_bar = theta_bar + del_theta;
    
    if norm(del_theta) < del_theta_tol
        break;
    end
end


function x_dot = state_equation(x, u)
    theta = x(1);
    theta_dot = x(2);
    
    g = 9.81;
    L = 2;
    m = 1;
    theta_ddot = -g/L*sin(theta) + 1/(m*L*L)*u;
    
    x_dot = [theta_dot;
             theta_ddot];
end

function z = observation_equation(x, m)
    offset = m(1);
    theta = x(1);
    
    z = offset + theta;
end

function [A, b] = build_A_b(theta, spline, t, z, R_z)
    [A_z, b_z] = build_A_b_z(theta, t, z, R_z, spline);
    [A_m, b_m] = build_A_b_m();
    [A_x, b_x] = build_A_b_x();
    [A_u, b_u] = build_A_b_u();
    
    A = sum(A_z) + A_m + A_x + A_u;
    b = sum(b_z) + b_m + b_x + b_u;
end

function [A_z, b_z] = build_A_b_z(theta, t, z, R_z, spline)
    % Now, there'll be matrices for each measurement, so these will be 3D
    A_z = zeros(length(theta), length(theta), size(z,2));
    b_z = zeros(length(theta), size(z,2));
    
    m = theta(end);
    
    for iii = 1:size(z,2)
        [J_x, J_m] = observation_jacobian(x, m);
        
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
    
    E_x = [phi_0, zeros(length(m))];
    e_x_bar = spline.evaluate(t_0, 0) - x_0_hat;
    
    A_x = (E_x')*inv(P_x)*E_x;
    b_x = (E_x')*inv(P_x)*e_x_bar;
end

function [A_u, b_u] = build_A_b_u(theta, spline, Q)
    
    m = theta(end);
    
    J = @(t) state_jacobian(x(t), u(t));
    
    E_u = @(t) [spline.get_phi(t, 1) - J(t)*spline.get_phi(t, 0), zeros(length(m))];
    e_u_bar = @(t) spline.evaluate(t, 1) - state_equation(spline.evaluate(t, 0), u(t));
    
    % Now we need to integrate two things to get A and B
    integrand_A = @(t) (E_u(t)')*inv(Q)*E_u(t);
    integrand_b = @(t) (E_u(t)')*inv(Q)*e_u_bar(t);
    
    ttt = t_0:0.01:t_f;
    AAA = zeros(length(theta));
    bbb = zeros(length(theta),1);
    
    for iii = 1:ttt
        AAA(:,:,iii) = integrand_A(ttt);
        bbb(:,iii) = integrand_b(ttt);
    end
    
    A_u = trapz(ttt, AAA, 3);
    b_u = trapz(ttt, bbb, 2);
end

function J = finite_difference_jacobian(f, x)
    del = 1e-5;
    
    f_0 = f(x);
    J = zeros(length(f_0), length(x));
    
    for iii = 1:length(x)
        del_x = zeros(size(x));
        del_x(iii) = del;
        del_f = f(x + del_x);
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

function J = objective_function(theta, t, z, R_z, spline, m_0, P_m, x_0, P_x, u, Q)
    c = theta(1:(end - 1));
    m = theta(end);
    
    spline.set_coefficients(c);
    
    e_z = zeros(2, length(t));
    
    for iii = 1:length(t)
        e_z(:,iii) = z(:,iii) - observation_equation(spline.evaluate(t(iii), 0), m);
    end
    
    J_z = zeros(1, length(t));
    
    for iii = 1:length(t)
        e_z_i = e_z(:,iii);
        J_z(iii) = 1/2*(e_z_i')*inv(R_z)*e_z_i;
    end
    
    e_m = m - m_0;
    J_m = 1/2*(e_m')*inv(P_m)*e_m;
    
    e_x = spline.evaluate(spline.a, 0) - x_0;
    J_x = 1/2*(e_x')*inv(P_x)*e_x;
    
    e_u = @(t) spline.evaluate(t, 1) - state_equation(spline.evaluate(t, 0), u(t));
    integrand = @(t) (e_u(t)')*inv(Q)*e_u(t);
    
    ttt = spline.a:0.1:spline.b;
    fff = zeros(size(ttt));
    for iii = 1:length(ttt)
        fff(iii) = integrand(ttt(iii));
    end
    
    J_u = 1/2*trapz(ttt, fff);
    
    J = J_m + J_x + sum(J_z) + J_u;
end
