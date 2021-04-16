function plot_wheel_robot(state_k)
clf;
hold on;

x = state_k(1,:);
y = state_k(2,:);
theta = state_k(3,:);

plot(x, y, '-r', 'LineWidth', 1);


ds = sqrt(mean(diff(x)).^2 + mean(diff(y)).^2);
num_downsample = floor(2/ds);

x_ds = downsample(x, num_downsample);
y_ds = downsample(y, num_downsample);
theta_ds = downsample(theta, num_downsample);

for iii = 1:length(x_ds)
    plot_robot(x_ds(iii), y_ds(iii), theta_ds(iii));
end

daspect([1,1,1]);
axis equal;
grid on;
end

function plot_robot(x, y, theta)
    % Draw the stick
    draw_rectangle(x, y, 0.1, 1.4, theta, [0, 0, 0]);
    % Draw the two wheels
    x_r = x + 1/2*cos(pi/2 - theta);
    y_r = y - 1/2*sin(pi/2 - theta);
    draw_rectangle(x_r, y_r, 0.75, 0.2, theta, [0, 0, 1]);
    x_l = x - 1/2*cos(pi/2 - theta);
    y_l = y + 1/2*sin(pi/2 - theta);
    draw_rectangle(x_l, y_l, 0.75, 0.2, theta, [0, 0, 1]);
end

function draw_rectangle(x, y, L, H, theta, rgb)
    R= ([cos(theta), -sin(theta); sin(theta), cos(theta)]);
    X=([-L/2, L/2, L/2, -L/2]);
    Y=([-H/2, -H/2, H/2, H/2]);
    T = R*[X; Y];
    x_lower_left=x+T(1,1);
    x_lower_right=x+T(1,2);
    x_upper_right=x+T(1,3);
    x_upper_left=x+T(1,4);
    y_lower_left=y+T(2,1);
    y_lower_right=y+T(2,2);
    y_upper_right=y+T(2,3);
    y_upper_left=y+T(2,4);
    x_coor=[x_lower_left x_lower_right x_upper_right x_upper_left];
    y_coor=[y_lower_left y_lower_right y_upper_right y_upper_left];
    patch('Vertices',[x_coor; y_coor]','Faces',[1 2 3 4],'Edgecolor','k','Facecolor',rgb,'Linewidth',1.2);
    axis equal;

end