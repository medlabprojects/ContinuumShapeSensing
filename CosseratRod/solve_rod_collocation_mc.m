function [p, q, p_cov, q_cov] = solve_rod_collocation_mc(f, f_cov, l, l_cov, K, N, M, L, p_0, R_0)

num_sims = 10;
dt = 0.1;

s_s = linspace(0, L, 1000);
pq_sims = zeros(7, length(s_s), num_sims);

figure(1);
clf;
hold on;
    
for iii = 1:num_sims
    f_process = sample_normal_process(f_cov, [0, L], dt);
    f_noise = @(t) f_process.evaluate(t);
    
    f_i = @(t) f(t) + f_noise(t);
    
    l_process = sample_normal_process(l_cov, [0, L], dt);
    l_noise = @(t) l_process.evaluate(t);
    
    l_i = @(t) l(t) + l_noise(t);
    
    [p_i, q_i] = solve_rod_collocation(f_i, l_i, K, N, M, L, p_0, R_0);
    pq_sims(1:3,:,iii) = p_i(s_s);
    pq_sims(4:7,:,iii) = q_i(s_s);
    
    if iii == 1
        plot_static_rod(p_i, q_i, f, l, L);
    else
        plot_static_rod(p_i, q_i, @(t) zeros(3, length(t)), @(t) zeros(3, length(t)), L);
        drawnow();
    end
end

% Now we have all of the simulated functions, lets fit the mean and cov
t_s = linspace(0, L, 1000);

p_1 = p_sims{1:3:end}


end