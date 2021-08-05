%% Use the shooting method to take some data
clear;

L = 1;
a = 1;
f = @(s) 25.*normpdf(s, 0.4, 0.025) + -10.*normpdf(s, 0.7, 0.025);

[s, x_truth] = solve_rod_shooting(a, L, f);

figure(1);
clf;

subplot(2,1,1);
plot(s, f(s));
title('Input Force');
xlabel('Arclength');
ylabel('Distrubuted Load');
subplot(2,1,2);
hold on;

plot(x_truth(:,1), x_truth(:,2));
title('Shape');
xlabel('x');
ylabel('y');

daspect([1,1,1]);

%% Subsample data to define y vector
ds = 100;

s_K = downsample(s, ds);

% Now we sample the x and y position to get y

y_tmp = downsample(x_truth, ds);
y_tmp = y_tmp(:,1:2);
y = reshape(y_tmp', [], 1);

plot(y(1:2:end), y(2:2:end), '.k');

%% Implementing continuous-arclength estimation

x_op = zeros(7*K,1);

% How do we initialize all of the wedges?
%   
% Should we do a linear interpolation for the Gamma in get_phi?
% Explain how I plan to implement the "measurement" of base and tip states.
%   Should I use if blocks for the g function on these?

max_iter = 100;

for iii = 1:max_iter
    
    phi_h
    % Return a function handle for the interpolation of x_op
    x_op_h = gp_interp(t, x_op, x_wed, phi_h, v_h, L_h, Q_c);
    
    F_h = @(t) df_dx(x_op_h(t));
    G_h = @(t) df_dw(x_op_h(t));
    v_h = @(t) f(x_op_h(t), u(t), 0);
    
    % Return a function handle for phi with two inputs t and t'
    phi_h = get_phi(F_h, t_0, t_f);
    
    % Build the required matrices for the linear solve
    v = build_v(x_wed_0, phi_h, v_h, t);
    F_inv = build_F_inv(phi_h, t);
    Q_inv = build_Q_inv(P_wed_0, phi_h, L_h, Q_c, t);
    
    R_inv = build_R_inv(); % This will be a blkdiag of all of the R's. They'll be different for different types of measurements.
    G = build_G(G_h, t);
    g = build_g(x_op_h, t);
    
    % Solve AAA*del_x = bbb
    AAA = (F_inv')*Q_inv*F_inv + (G')*R_inv*G;
    bbb = (F_inv')*Q_inv*(v - F_inv*x_op) + (G')*R_inv*(y - g);
    
    del_x = AAA \ bbb;
    
    x_op = x_op + del_x;
    
    if norm(x_op) < 1e-6
        x_hat = x_op;
        break;
    end
end

% Now use GP interpolation to compute the estimate at other times of
% interest, x_hat_tau, P_hat_tautau

function x_dot = f_h(x)
    p_x = x(1);
    p_y = x(2);
    theta = x(3);
    m = x(4);
    n_x = x(5);
    n_y = x(6);
    force = x(7);
    
    p_x_dot = cos(theta);
    p_y_dot = sin(theta);
    theta_dot = -a*m;
    m_dot = sin(theta)*n_x + cos(theta)*n_y;
    n_x_dot = -sin(theta)*force;
    n_y_dot = -cos(theta)*force;
    force_dot = 0;
    
    x_dot = [p_x_dot; p_y_dot; theta_dot; m_dot; n_x_dot; n_y_dot; force_dot];
end

function F = df_dx(x)
    p_x = x(1);
    p_y = x(2);
    theta = x(3);
    m = x(4);
    n_x = x(5);
    n_y = x(6);
    force = x(7);
    
    F = [0, 0, -sin(theta),                          0, 0,          0,          0
         0, 0,  cos(theta),                          0, 0,          0,          0
         0, 0,  0,                                  -a, 0,          0,          0
         0, 0,  n_x*cos(theta) - n_y*sin(theta),     0, sin(theta), cos(theta), 0
         0, 0, -force*cos(theta),                    0, 0,          0,         -sin(theta)
         0, 0,  force*(sin(theta)),                  0, 0,          0,         -cos(theta)
         0, 0,  0,                                   0, 0,          0,          0];
end

function L = df_dw()
    L = [0; 0; 0; 0; 0; 0; 1];
end

function get_phi(F)
    



end

function gp_interp(x, tau)



end

