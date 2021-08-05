clear;

% Time interval a to b
a = 0;
b = 6;

num_segments = 10;
segment_length = (b - a)/num_segments;

num_knots = num_segments + 7;
num_basis = num_segments + 3;

% The dimension of the space the spline maps to
output_dim = 3;

knot_sequence = (a - 3*segment_length):segment_length:(b + 3*segment_length);
coefficients = zeros(output_dim*num_basis,1);

coefficients(13) = 1;
coefficients(14) = 1;
coefficients(15) = 1;

% coefficients = rand(size(coefficients));

t = linspace(knot_sequence(1), knot_sequence(end), 1000);
x = zeros(3, length(t));
x_dot = zeros(3, length(t));
x_ddot = zeros(3, length(t));

for iii = 1:length(t)
    if (t(iii) <= a) || (t(iii) >= b)
        x(:,iii) = nan(3,1);
        x_dot(:,iii) = nan(3,1);
        x_ddot(:,iii) = nan(3,1);
    else
        phi = get_phi(t(iii), 0, knot_sequence, num_basis, output_dim);
        phi_dot = get_phi(t(iii), 1, knot_sequence, num_basis, output_dim);
        phi_ddot = get_phi(t(iii), 2, knot_sequence, num_basis, output_dim);
        
        x(:,iii) = phi*coefficients;
        x_dot(:,iii) = phi_dot*coefficients;
        x_ddot(:,iii) = phi_ddot*coefficients;
    end
end

figure(1);
clf;
hold on;
grid on;

plot(t, x);

title('Position vs. Time');

figure(2);
clf;
hold on;
grid on;

plot(t, x_dot);
plot(t(2:end), diff(x,1,2)./mean(diff(t)));

title('Velocity vs. Time with Finite Diff');

figure(3);
clf;
hold on;
grid on;

plot(t, x_ddot);
plot(t(3:end), diff(x,2,2)./(mean(diff(t))^2));

title('Accelleration vs. Time with Finite Diff');

function phi = get_phi(t, deriv_order, knot_sequence, num_basis, output_dim)
    
    B = 1/6*[1  4  1  0
            -3  0  3  0
             3 -6  3  0
            -1  3 -3  1];
    
    D = [0 1 0 0
         0 0 2 0
         0 0 0 3
         0 0 0 0];
    
    phi = zeros(output_dim, num_basis*output_dim);
    interval = find_interval(t, knot_sequence);
    
    t_ii = knot_sequence(interval);
    t_ip1 = knot_sequence(interval + 1);
    u = (t - t_ii)/(t_ip1 - t_ii);
    U = [1, u, u^2, u^3];
    
    for iii = 1:output_dim
        b1 = output_dim*(interval - 4) + iii;
        b2 = output_dim*(interval - 3) + iii;
        b3 = output_dim*(interval - 2) + iii;
        b4 = output_dim*(interval - 1) + iii;
        
        P_i = zeros(4, num_basis*output_dim);
        P_i(1, b1) = 1;
        P_i(2, b2) = 1;
        P_i(3, b3) = 1;
        P_i(4, b4) = 1;
        
        phi(iii,:) = U*(D^deriv_order)*B*P_i;
    end
    
    segment_length = mean(diff(knot_sequence));
    phi = phi/(segment_length^deriv_order);
end

function interval = find_interval(t, t_i)
    interval = find(diff(t >= t_i));
end