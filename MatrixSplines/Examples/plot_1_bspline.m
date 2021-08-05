clear;

B = 1/6*[1  4  1  0
        -3  0  3  0
         3 -6  3  0
        -1  3 -3  1];
    
t_i = -3:9;
c = zeros(9,1);
c(4) = 1;

%% For the interval between 0 and 1, we're using c-3, c-2, c-1, c0
t01 = linspace(0, 1, 100);
x01 = zeros(size(t01));
c01 = c(1:4);

for iii = 1:length(t01)
    u = (t01(iii) - 0)/(1 - 0);
    U = [1, u, u^2, u^3];
    x01(iii) = U*B*c01;
end

figure(1);
clf;
hold on;

plot(t01, x01);

%% For the interval between 1 and 2, we're using c-2, c-1, c0, c1
t12 = linspace(1, 2, 100);
x12 = zeros(size(t12));
c12 = c(2:5);

for iii = 1:length(t01)
    u = (t12(iii) - 1)/(2 - 1);
    U = [1, u, u^2, u^3];
    x12(iii) = U*B*c12;
end

plot(t12, x12);

%% For the interval between 2 and 3, we're using c-1, c0, c1, c2
t23 = linspace(2, 3, 100);
x23 = zeros(size(t23));
c23 = c(3:6);

for iii = 1:length(t01)
    u = (t23(iii) - 2)/(3 - 2);
    U = [1, u, u^2, u^3];
    x23(iii) = U*B*c23;
end

plot(t23, x23);

% For the interval between 3 and 4, we're using c0, c1, c2, c3
t34 = linspace(3, 4, 100);
x34 = zeros(size(t34));
c34 = c(4:7);

for iii = 1:length(t01)
    u = (t34(iii) - 3)/(4 - 3);
    U = [1, u, u^2, u^3];
    x34(iii) = U*B*c34;
end

plot(t34, x34);





function interval = find_interval(t, t_i)
    t > t_i
end