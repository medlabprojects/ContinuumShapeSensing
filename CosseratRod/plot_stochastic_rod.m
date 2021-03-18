function plot_stochastic_rod(s, state, f)

% First convert forces from local frame to global
fLocal = f(s)';
fGlobal = localForceToGlobal(fLocal, state);

x = state(:,1);
y = state(:,2);
z = state(:,3);

figure(1);
clf;
hold on;

plot3(x, y, z, '-k', 'LineWidth', 5);
quiver3(x, y, z, fGlobal(:,1), fGlobal(:,2), fGlobal(:,3));

axis tight;
grid on;
view(-120,30);
daspect([1 1 1]);

xlabel('x')
ylabel('y')
zlabel('z')
end
    
function fGlobal = localForceToGlobal(fLocal, state)
    N = size(state, 1);
    fLocal = [fLocal, zeros(N, 1)];
    fGlobal = zeros(size(fLocal));
    
    for iii = 1:N
        [~, R] = UnpackRodState(state(iii,:)');
        fGlobal(iii,:) = (R*(fLocal(iii,:)'))';
    end
end
