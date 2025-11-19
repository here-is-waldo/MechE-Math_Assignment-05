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


    % Uncomment the following code to plot:

    
    % Vx_linear = V_lin(:,1)-V_eq(1);
    % Vy_linear = V_lin(:,2)-V_eq(2);
    % Vtheta_linear = V_lin(:,3)-V_eq(3);
    % 
    % figure(2);
    % subplot(2,1,1)
    % title('Linearized vs Nonlinear Dynamics, e = 10^-^3');
    % hold on
    % 
    % plot(t_nl, Vx_linear,'k','linewidth',2)
    % plot(t_nl, Vy_linear,'k','linewidth',2)
    % plot(t_nl, Vtheta_linear,'k','linewidth',2)
    % plot(t_nl, V_nl(:,1)-V_eq(1),'r--','linewidth',2)
    % plot(t_nl, V_nl(:,2)-V_eq(2),'b--','linewidth',2)
    % plot(t_nl, V_nl(:,3)-V_eq(3),'g--','linewidth',2)
    % legend("X (Linear)", "Y (Linear)", "theta (Linear)", "X (Nonlinear)", "Y (Nonlinear)", "theta (Nonlinear)")
    % xlabel('Time (s)'); ylabel('Displacement from Equilibrium');
    % legend('Location','northeast');
    % grid on;
    % 
    % 
    % % hold off
    % 
    % subplot(2,1,2)
    % hold on
    % title('Linearized vs Nonlinear Dynamics, e = 1');
    % epsilon = 1;
    % V0 = V_eq + epsilon * [1; 0; 0; 0; 0; 0];
    % 
    % tspan = [0 5];
    % opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
    % 
    % [t_lin, V_lin] = ode45(linear_rate, tspan, V0, opts);
    % [t_nl,  V_nl]  = ode45(nonlinear_rate, tspan, V0, opts);
    % 
    % Vx_linear = V_lin(:,1)-V_eq(1);
    % Vy_linear = V_lin(:,2)-V_eq(2);
    % Vtheta_linear = V_lin(:,3)-V_eq(3);
    % 
    % figure(2);
    % subplot(2,1,2)
    % title('Linearized vs Nonlinear Dynamics, e = 1');
    % hold on
    % 
    % plot(t_lin, Vx_linear,'k','linewidth',2)
    % hold on
    % plot(t_lin, Vy_linear,'k','linewidth',2)
    % plot(t_lin, Vtheta_linear,'k','linewidth',2)
    % plot(t_nl, V_nl(:,1)-V_eq(1),'r--','linewidth',2)
    % plot(t_nl, V_nl(:,2)-V_eq(2),'b--','linewidth',2)
    % plot(t_nl, V_nl(:,3)-V_eq(3),'g--','linewidth',2)
    % legend("X (Linear)", "Y (Linear)", "theta (Linear)", "X (Nonlinear)", "Y (Nonlinear)", "theta (Nonlinear)")
    % xlabel('Time (s)'); ylabel('Displacement from Equilibrium');
    % legend('Location','northeast');
    % grid on;


end
