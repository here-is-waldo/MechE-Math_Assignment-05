function V_eq = find_equilibrium(box_params, V_guess)

    
    rate_func = @(V) box_rate_func(0, V, box_params);
    fprintf('Finding equilibrium configuration...\n');
    V_eq = f_multi_newton(rate_func, V_guess);
    fprintf('Equilibrium state found:\n');
    disp(V_eq.');
    fprintf('Verifying equilibrium with ODE45...\n');
    tspan = [0 5];
    my_rate_func = @(t,V) box_rate_func(t,V,box_params);

    opts = odeset('RelTol',1e-8,'AbsTol',1e-9);
    [t_list, V_list] = ode45(my_rate_func, tspan, V_eq, opts);

    figure;
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
