function plot_static_rod(x, f, l, L)

% First convert forces from local frame to global
s = linspace(0, L, 1000);
xs = x(s);
fs = f(s);

f_global = local_force_to_global(fs, xs);

x = xs(1,:);
y = xs(2,:);
z = xs(3,:);

f_x = f_global(1,:);
f_y = f_global(2,:);
f_z = f_global(3,:);

hold on;

plot3(x, y, z, '-k', 'LineWidth', 5);
quiver3(x, y, z, f_x, f_y, f_z);

axis tight;
grid on;
view(-120,30);
daspect([1 1 1]);

xlabel('x')
ylabel('y')
zlabel('z')
end
    
function f_global = local_force_to_global(f_local, x)
    N = size(x, 2);
    f_global = zeros(size(f_local));
    
    for iii = 1:N
        [~, q] = unpack_rod_state(x(:,iii));
        R = quat2rotm(q');
        f_global(:,iii) = R*f_local(:,iii);
    end
end