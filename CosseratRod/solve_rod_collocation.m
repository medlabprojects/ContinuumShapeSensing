function x = solve_rod_collocation(f, l, K, N, M, L, p_0, R_0)
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
x_spline = BSpline([0, L], d, 13, k);
x_spline.C(4,:) = 1; % For quaternions
C0 = x_spline.C;

s = linspace(0, L, N);

options = optimoptions('fsolve');
options.Display = 'iter';
options.Algorithm = 'levenberg-marquardt';

obj = @(C) compute_total_residual(C, x_spline, f, l, s, K, p_0, R_0);
C = fsolve(obj, C0, options);

x_spline.C = C;

x = @(s) x_spline.evaluate(s);
end

function res = compute_total_residual(C, x_spline, f, l, s, K, p_0, R_0)
    % Update spline and get derivative too
    x_spline.C = C;
    x_spline_deriv = x_spline.deriv();
    
    % Evaluate on the collocation points
    x_s = x_spline.evaluate(s);
    xp_s = x_spline_deriv.evaluate(s);
    
    f_s = f(s);
    l_s = l(s);
    
    % Determine the goodness of fit for the differential eq. constraint
    xp_res = zeros(13, length(s));
    
    e_3 = [0; 0; 1];
    
    for iii = 1:length(s)
        [~, qi, ni, mi] = unpack_rod_state(x_s(:,iii));
        fi = f_s(:,iii);
        li = l_s(:,iii);
        
        Ri = quat2rotm(qi');
        ui = (K\(Ri'))*mi;
        
        p_dot = Ri*e_3;
        q_dot = 1/2*quatmultiply(qi', [0; ui]');
        n_dot = -fi;
        m_dot = -cross(p_dot, ni) - li;
        
        xp_i = [p_dot; q_dot'; n_dot; m_dot];
        xp_res(:,iii) = xp_i - xp_s(:,iii);
    end
    
    % Determine goodness of fit for the base and tip constraints
    x_0 = x_s(:,1);
    [p_0_coll, q_0_coll] = unpack_rod_state(x_0);
    
    base_res = [p_0_coll - p_0; [1; 0; 0; 0] - quatmultiply(q_0_coll', rotm2quat(R_0'))'];
    
    x_L = x_s(:,end);
    tip_res = x_L(8:end);
    
    % Determine the goodness of fit for the quaternions
    q = x_s(4:7,:);
    q_res = quatnormalize(q')' - q;

    res = [xp_res(:); base_res(:); tip_res(:); q_res(:)]; 
end