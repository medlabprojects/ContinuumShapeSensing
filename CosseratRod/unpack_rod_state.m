function [p, R, n, m] = unpack_rod_state(x)
p = x(1:3,1);
R = [x(4:6), x(7:9), x(10:12)];
n = x(13:15);
m = x(16:18);

R = orthogonalize(R);
end