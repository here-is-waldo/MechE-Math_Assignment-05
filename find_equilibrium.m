%% locate an equilibrium configuration

% INPUTS:
%   box_params: a struct containing the parameters that describe the system
%   Fields:
%       box_params.m: mass of the box
%       box_params.I: moment of inertia w/respect to centroid
%       box_params.g: acceleration due to gravity
%       box_params.k_list: list of spring stiffnesses
%       box_params.l0_list: list of spring natural lengths
%       box_params.P_world: 2 x n list of static mounting
%                       points for the spring (in the world frame)
%       box_params.P_box: 2 x n list of mounting points
%                       for the spring (in the box frame)
%   V_guess = [x ; y; theta; dxdt; dydt; dthetadt]: state vector eq. guess
% OUTPUTS
%   V_eq = [x ; y; theta; dxdt; dydt; dthetadt]: state vector for the
%                       equlibirum position

function V_eq = find_equilibrium(box_params, V_guess)

    
    rate_func = @(V) box_rate_func(0, V, box_params);
    fprintf('Finding equilibrium configuration...\n');
    solver_params = struct();
    V_eq = multi_newton(rate_func, V_guess, solver_params);
    fprintf('Equilibrium state found:\n');
    disp(V_eq.');
    fprintf('Verifying equilibrium with ODE45...\n');
    tspan = [0 5];
    my_rate_func = @(t,V) box_rate_func(t,V,box_params);

    opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
    [t_list, V_list] = ode45(my_rate_func, tspan, V_eq, opts);


    figure(1);
    plot(t_list, V_list(:,1)-V_eq(1), 'r', 'LineWidth',1.2);
    hold on;
    plot(t_list, V_list(:,2)-V_eq(2), 'b', 'LineWidth',1.2);
    plot(t_list, V_list(:,3)-V_eq(3), 'k', 'LineWidth',1.2);
    xlabel('Time (s)');
    ylabel('Deviation from Equilibrium');
    legend('\Delta x','\Delta y','\Delta \theta');
    title('Equilibrium Verification (ODE45)');
    grid on;
    

    fprintf('If the system remains stationary (flat lines), equilibrium is valid.\n');
end
