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

