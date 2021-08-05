function x_dot = state_deriv(s, x, f, a)
p_x = x(1);
p_y = x(2);
theta = x(3);
m = x(4);
n_x = x(5);
n_y = x(6);

p_x_dot = cos(theta);
p_y_dot = sin(theta);

theta_dot = -a*m;

m_dot = sin(theta)*n_x + cos(theta)*n_y;

n_x_dot = -sin(theta)*f(s);
n_y_dot = -cos(theta)*f(s);

x_dot = [p_x_dot; p_y_dot; theta_dot; m_dot; n_x_dot; n_y_dot];

end