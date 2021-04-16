% Test Two Wheeled Robot Example

T = 10; % seconds
dt = 0.01;

t_k = 0:dt:T;

% Now we're going to assume some inputs
v_r = 1 + 0.25*t_k + 0.5*sin(t_k);
v_l = 1 + 0.25*t_k + 0.5*cos(t_k);

v_k = [v_l; v_r];
state_0 = [1; 1; deg2rad(45)];

state_k = sim_wheel_robot(state_0, t_k, v_k);

figure(1);
plot_wheel_robot(state_k);
