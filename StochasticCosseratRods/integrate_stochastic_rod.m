function [s, x, x_cov] = integrate_stochastic_rod(x_0, x_cov_0, f, f_cov, l, l_cov, N, s_span, K)
% Given an intialization for state x_0 and associated covariance matrix 
% x_cov_0, the number of points in the integration N, the initial and final 
% arclengths s_span, and arclength parameterized mean and covariance 
% distribited load function handles (f, f_cov, l, l_cov), and the rod 
% bending matrix K, this returns the mean x at every integration point s.
% Furthermore, the state covariance x_cov at every s is also returned.

% The rod state x is defined by:
%   x_1-3:   position of rod, p
%   x_4-12:  rotation matrix of rod, R
%   x_13-15: internal force of rod, n
%   x_16-18: internal moment of rod, m

process = @(s, y) kalman_bucy(s, y, f, f_cov, l, l_cov, K, Q);
y_0 = pack_y(x_0, x_cov_0);

[s, y]  = ode45(process, linspace(s_span(1), s_span(2), N), y_0);

[x, x_cov] = unpack_y(y);
end

function y_dot = kalman_bucy(s, y, f, f_cov, l, l_cov, K, Q)
    [x, x_cov] = unpack_y(y);
    
    % Compute the mean of x_dot
    x_dot = simplified_cosserat_model(s, x, f, l, K);
    
    % Evaluate distributed load functions at s
    f_s = f(s);
    f_cov_s = f_cov(s);
    l_s = l(s);
    l_cov_s = l_cov(s);
    
    % Compute Jacobians for uncertainty propogation
    F_x = compute_jacobian(@(x) simplified_cosserat_model(s, x, f, l, K), x, 1e-5);
    F_f = compute_jacobian(@(f) simplified_cosserat_model(s, x, f, l, K), f_s, 1e-5);
    F_l = compute_jacobian(@(l) simplified_cosserat_model(s, x, f, l, K), l_s, 1e-5);
    
    % Known GP parameters Q
    x_cov_dot = Q;
    
    % Propogate state uncertainty
    x_cov_dot = x_cov_dot + F_x*x_cov + x_cov*F_x;
    
    % Propogate force uncertainty
    x_cov_dot = x_cov_dot + F_f*f_cov_s + f_cov_s*F_f;
    
    % Propogate moment uncertainty
    x_cov_dot = x_cov_dot + F_l*l_cov_s + l_cov_s*F_l;
    
    y_dot = pack_y(x_dot, x_cov_dot);
end

function [x, x_cov] = unpack_y(y)
    length_x = 18;
    x = y(1:length_x);
    x_cov_vec = y((length_x + 1):end);
    x_cov = reshape(x_cov_vec, length_x, length_x);
end

function y = pack_y(x, x_cov)
    y = [x; x_cov(:)];
end

function x_dot = simplified_cosserat_model(s, x, f, l, K)
    [~, R, n, m] = unpack_rod_state(x);
    
    e_3 = [0; 0; 1];
    u = (K\(R'))*m;
    
    p_dot = R*e_3;
    R_dot = R*hat(u);
    
    n_dot = -f(s);
    m_dot = -cross(p_dot, n) - l(s);
    
    x_dot = pack_rod_state(p_dot, R_dot, n_dot, m_dot);
end

function matrix = hat(vector)
    x = vector(1);
    y = vector(2);
    z = vector(3);
    
    matrix = [0 -z y; z 0 -x; -y x 0];
end

function J = compute_jacobian(f, x, del)
    f0 = f(x);
    J = zeros(size(f0), size(x));
    
    for iii = 1:length(x)
        x_left = x;
        x_right = x;
        
        x_left(iii) = x_left(iii) - del;
        x_right(iii) = x_right(iii) + del;
        
        df = f(x_right) - f(x_left);
        dx = 2*del;
        
        J(:,iii) = df./dx;
    end
end

