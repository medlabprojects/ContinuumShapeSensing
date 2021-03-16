function [s, state] = solve_stochastic_rod(p_0, p_cov_0, R_0, R_cov_0, D, E, G, N, sSpan, f, f_cov, l, l_cov)
% We assume here that both rod ends are force and moment free. We will
% initialize all states and forcing functions and then solve the boundary
% value problem subject to the boundary conditions using an optimization
% method.

% Problem unknowns: initial covariances p_0_cov, R_0_cov and posterior
% loading distributions f posterior, l posterior

% Problem knowns: internal forces and moments at the endpoints MUST be
% certain. We are considering the case of only distributed loads which
% should in principle encompass discrete loading too. Therefore, 
%   n_0 = m_0 = 0, n_0_cov = m_0_cov = 0,
%   n_L = m_L = 0, n_L_cov = n_L_cov = 0.

% Check: the above may be inconsistant. I.e., we may be created an
% ill-posed problem and should instead not consider all the covariances as
% unknown? We may need to look into collocation methods with B-splines for 
% this problem. I think it can be recast... Maybe the forward and inverse
% problems can both be recast as regularized optimization problems?

% Solution approach: make these known conditions as true as possible by 
% varying the unknown functions and parameters. The unknown functions are 
% parameterized by B-spline coefficients. All unknowns (including function
% parameters) are grouped into a vector for optimization. Constraints are 
% setup for the known endpoint conditions. The cost is the deviation of the
% loading functions from nominal inner product-ed with their gaussian 
% process gains. With constraints, cost, and variables all defined, just 
% use fsolve 

All of the parameters are varied 
% to produce the most likely feasible solution in continuous arclength.

p_0
R_0

x_0 = pack_rod_state(p_0, R_0, n_0, m_0)









I = pi/4*(D/2)^4; %m^4
K = [E*I 0 0; 0 E*I 0; 0 0 2*G*I];

n_0 = zeros(3,1);
u_0 = zeros(3,1);

state_0 = PackRodState(p_0, R_0, n_0, u_0);
nu_0_guess = zeros(6,1);

options = optimoptions('fsolve');
options.MaxFunctionEvaluations = 1e4;

nu_0 = fsolve(@(nu_0) compute_rod_tip_residual(nu_0, state_0, N, sSpan, f, K), nu_0_guess, options);

state_0 = PackRodState(p_0, R_0, nu_0(1:3), nu_0(4:6));

[s, state] = IntegrateStaticRod(state_0, N, sSpan, f, K);
end


function res = compute_rod_tip_residual(nu_0, state_0, N, sSpan, f, K)
    % Unknown is the internal forces and moments at rod base
    state_0(13:18) = nu_0;
    
    [~, state] = IntegrateStaticRod(state_0, N, sSpan, f, K);
    
    % Want the internal forces and moments to be zero at rod tip
    [~, ~, n_f, u_f] = UnpackRodState(state(end,:)');
    res = [n_f; u_f];
end