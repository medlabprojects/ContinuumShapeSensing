function [s, x] = solve_rod_shooting(a, L, f)

obj = @(y) end_condition(y, f, a, L);
y_0 = [0;0;0];

y = fsolve(obj, y_0);

[c, s, x] = obj(y);

if norm(c) > 0.01
    warning('Shooting method did not converge.');
end

end


function [c, s, x] = end_condition(y, f, a, L)

m_0 = y(1);
n_x_0 = y(2);
n_y_0 = y(3);

s = linspace(0, L, 1000);

x_0 = [0; 0; 0; m_0; n_x_0; n_y_0];

deriv = @(s, x) state_deriv(s, x, f, a);

[s, x] = ode45(deriv, s, x_0);

c = x(end,4:6);
end