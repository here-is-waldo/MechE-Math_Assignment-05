function [A, J_approx] = linearize_system(box_params, V_eq)

    J_approx = approximate_jacobian(@(V) box_rate_func(0, V, box_params), V_eq);
    A = J_approx;

    linear_rate = @(t, V) A * (V - V_eq);
    nonlinear_rate = @(t, V) box_rate_func(t, V, box_params);

    epsilon = 1e-3;
    V0 = V_eq + epsilon * [1; 0; 0; 0; 0; 0];

    tspan = [0 5];
    opts = odeset('RelTol',1e-8,'AbsTol',1e-9);

    [t_lin, V_lin] = ode45(linear_rate, tspan, V0, opts);
    [t_nl,  V_nl]  = ode45(nonlinear_rate, tspan, V0, opts);

    figure(); hold on;
    plot(t_nl, V_nl(:,1)-V_eq(1), 'r', 'LineWidth',1.2, ...
         'DisplayName','Nonlinear ∆x');
    plot(t_lin, V_lin(:,1)-V_eq(1), 'b--', 'LineWidth',1.2, ...
         'DisplayName','Linearized ∆x');
    xlabel('Time (s)'); ylabel('Displacement ∆x');
    legend('Location','best');
    grid on;
    title('Linearized vs Nonlinear Dynamics (ODE45 Verification)');
end
