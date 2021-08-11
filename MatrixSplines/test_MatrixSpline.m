clear;

%% First just do some quick position, velocity, accelleration tests

a = 0;
b = 5;
num_segments = 10;
output_dim = 3;

spline = MatrixSpline();
spline.setup(a, b, num_segments, output_dim);

coefficients = zeros(size(spline.coefficients));

coefficients(1) = -20;
coefficients(4) = 1;
coefficients(7) = 1;

% coefficients = rand(size(coefficients));
spline.set_coefficients(coefficients);

t = linspace(a, b, 1000);

x = spline.evaluate(t, 0);
x_dot = spline.evaluate(t, 1);
x_ddot = spline.evaluate(t, 2);

figure(1);
clf;
hold on;

plot(t, x);

title('Position vs. Time');

figure(2);
clf;
hold on;

plot(t, x_dot);
plot(t(2:end), diff(x,1,2)./mean(diff(t)));

title('Velocity vs. Time with Finite Diff');

figure(3);
clf;
hold on;

plot(t, x_ddot);
plot(t(3:end), diff(x,2,2)./(mean(diff(t))^2));

title('Accelleration vs. Time with Finite Diff');

%% Now we can do some interpolation

t_data = 0:1:5;
x_data = [sin(t_data); cos(t_data)];

spline = MatrixSpline();
spline.fit(t_data, x_data, 1);

t = linspace(spline.a, spline.b, 1000);
x = spline.evaluate(t, 0);
x_gt = [sin(t); cos(t)];

figure(4);
clf;
hold on;

plot(t_data, x_data, '.r', 'MarkerSize', 20);
plot(t, x_gt);
plot(t, x);

legend('','','','','');