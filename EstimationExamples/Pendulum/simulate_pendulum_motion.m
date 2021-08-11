function [t, x, z] = simulate_pendulum_motion(x_0, u, t_span, T, m, Q)

t = t_span(1):T:t_span(2);
x = zeros(2, length(t));
z = zeros(1, length(t));

x(:,1) = x_0;
z(:,1) = observation_equation(x_0, m);

for iii = 2:length(t)
    x_dot = state_equation(x(:,iii - 1), u(t(iii))) + mvnrnd([0;0], Q, 1)';
    x(:,iii) = x(:,iii - 1) + x_dot*T;
    z(:,iii) = observation_equation(x(:,iii), m);
end

end