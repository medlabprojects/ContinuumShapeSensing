function state_k = sim_wheel_robot(state_0, t_k, v_k)
    state_k = zeros(3, length(t_k));
    state_k(:,1) = state_0;
    
    dt = mean(diff(t_k));
    
    for iii = 2:length(t_k)
        state_dot = wheel_robot_dynamics(state_k(:,iii - 1), v_k(:,iii - 1));
        state_k(:,iii) = state_k(:,iii - 1) + state_dot*dt;
    end
end