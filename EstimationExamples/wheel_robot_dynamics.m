function state_dot = wheel_robot_dynamics(state, v)
v_l = v(1);
v_r = v(2);

v = 1/2*(v_l + v_r);

x = state(1);
y = state(2);
theta = state(3);

x_dot = cos(theta)*v;
y_dot = sin(theta)*v;
theta_dot = v_r - v_l;

state_dot = [x_dot; y_dot; theta_dot];
end