clear;

B = 1/6*[1  4  1  0
        -3  0  3  0
         3 -6  3  0
        -1  3 -3  1];

t_i = linspace(-3,9,13);

c = zeros(9,1);
c(4) = 1;
c(5) = -.7;
c(6) = 2;

t = linspace(0, 5, 1000);
x = zeros(size(t));

for iii = 1:length(t)
    interval = find_interval(t(iii), t_i);
    c_i = c((interval - 3):interval);
    t_ii = t_i(interval);
    t_ip1 = t_i(interval + 1);
    u = (t(iii) - t_ii)/(t_ip1 - t_ii);
    U = [1, u, u^2, u^3];
    x(iii) = U*B*c_i;
end

figure(1);
clf;
hold on;

plot(t, x);


function interval = find_interval(t, t_i)
    interval = find(diff(t >= t_i));
end