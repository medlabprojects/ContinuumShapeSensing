function [p, q, n, m] = unpack_rod_state(x)
p = x(1:3);
q = x(4:7);
n = x(8:10);
m = x(11:13);

% Unitize q
q = q./(norm(q));
end