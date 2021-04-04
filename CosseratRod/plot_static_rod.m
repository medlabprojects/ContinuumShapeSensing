function plot_static_rod(p, q, f, l, L)

% First convert forces from local frame to global
s = linspace(0, L, 1000);
p_s = p(s);
q_s = q(s);
f_s = f(s);

f_s = local_force_to_global(f_s, q_s);

x = p_s(1,:);
y = p_s(2,:);
z = p_s(3,:);

f_x = f_s(1,:);
f_y = f_s(2,:);
f_z = f_s(3,:);

hold on;

N = 100;
plot_transforms(downsample(p_s', N), downsample(q_s', N), L/30);

plot3(x, y, z, '-k', 'LineWidth', 5);
quiver3(x, y, z, f_x, f_y, f_z, 3.*L);

axis tight;
grid on;
view(-120,30);
daspect([1 1 1]);

xlabel('x')
ylabel('y')
zlabel('z')
end
    
function f_global = local_force_to_global(f_local, q)
    f_global = zeros(size(f_local));
    
    for iii = 1:size(q,2)
        R = quat2rotm(q(:,iii)');
        f_global(:,iii) = -(R')*f_local(:,iii);
    end
end

function plot_transforms(p, q, sz)
    R = quat2rotm(q);
    
    for iii = 1:size(p,1)
        draw_coordinate_system(sz, R(:,:,iii)', p(iii,:), 'rgb');
    end
end
