function [p, q] = solve_rod_collocation(f, l, K, N, M, L, p_0, R_0)
% Given the loading on the rod f and l (function handles), the rod diameter
% D, the elastic modulus E, and shear modulus G, the number of collocation
% points N, the length of the rod L, and initial, this function uses a 
% collocation method with B-splines to attempt to solve an approximating 
% solution for the state x (function handle). 

% We assume here that only the tip of the rod is force and moment free. The
% base of the rod can have any loading condition.

% We will initialize all states and forcing functions and then solve the 
% boundary value problem subject to the boundary conditions using an 
% optimization method.

% N is number of collocation points
% M is dimension of spline space which is different.

d = 5;
k = M - d - 1;

% Local curvature
u_spline = BSpline([0, L], d, 3, k);

% We're going to solve for the correct C
C0 = u_spline.C;

s = linspace(0, L, N);

options = optimoptions('fsolve');
options.Display = 'iter';
options.Algorithm = 'levenberg-marquardt';

obj = @(C) compute_residual(C, u_spline, f, l, s, K, p_0, R_0);
C = fsolve(obj, C0, options);

[~, p, q] = obj(C);
end

function [res, p, q, n, m] = compute_residual(C, u_spline, f, l, s, K, p_0, R_0)

    e_3 = [0; 0; 1];
    
    u_spline.C = C;
    u = @(s) u_spline.evaluate(s);
    v = @(s) e_3;
    
    % Forward integrate twist functions to get geometry information
    x_0 = [p_0; rotm2quat(R_0)'; 0];
    
    obj = @(s, x) deriv_pq(s, x, u, v);
    [~, x] = ode45(obj, s, x_0);
    x = x';
    
    if any(abs(x(end,:)) > 1e-3)
        warning('Quaternion residual greater than 0.001 during integration');
    end
  
    p_spline = BSpline.interpolate(x(1:3,:), s, 5);
    q_spline = BSpline.interpolate(x(4:7,:), s, 5);  
    p_dot_spline = p_spline.deriv();
    
    p = @(s) p_spline.evaluate(s);
    q = @(s) q_spline.evaluate(s);
    p_dot = @(s) p_dot_spline.evaluate(s);
    
    % Backward integrate to obtain mechanicals information
    y_0 = zeros(6,1);
    
    obj = @(s, y) deriv_nm(s, y, f, l, p_dot);
    [~, y] = ode45(obj, fliplr(s), y_0);
    y = fliplr(y');
    
    n_spline = BSpline.interpolate(y(1:3,:), s, 5);
    m_spline = BSpline.interpolate(y(4:6,:), s, 5);
    
    n = @(s) n_spline.evaluate(s);
    m = @(s) m_spline.evaluate(s);
    
    % Now, finally, we can close the loop on u and v
    res = u(s) - inv(K)*(quatrotate(q(s)', m(s)')');
%     v_res = v(s) - inv
end

function x_dot = deriv_pq(s, x, u, v)
    p = x(1:3);
    q = x(4:7);
    
    p_dot = quatrotate(q', v(s)')';
    q_dot = 1/2*quatmultiply(q', [0; u(s)]')';
    
    % Quaternion constraint = 0
    c = (q')*q - 1;
    
    x_dot = [p_dot; q_dot; c];
end


function y_dot = deriv_nm(s, y, f, l, p_dot)
    n = y(1:3);
%     m = y(4:6);
    
    n_dot = -f(s);
    m_dot = -cross(p_dot(s), n) - l(s);
    
    y_dot = [n_dot; m_dot];
end




