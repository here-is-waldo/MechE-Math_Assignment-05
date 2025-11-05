

function simulate_box()
    
    % get params
    box_params = get_Orion_params();
    
    % load the system parameters into the rate function
    % via an anonymous function
    my_rate_func = @(t_in, V_in) box_rate_func(t_in, V_in, box_params);
    
    % more initial conditions
    x0 = 2; y0 = 3; theta0 = 0.331;
    vx0 = 0.1; vy0 = 0.05; vtheta0 = 0.2;
    V0 = [x0; y0; theta0; vx0; vy0; vtheta0];
    tspan = [0, 5];

    % run the integration (using ode45 cuz my integrator isn't reliable)
    [tlist, Vlist] = ode45(my_rate_func, tspan, V0);

    % do the recording and animate stuff

    % initialize spring plot
    num_zigs = 6; w = 3;
    spring_plot_struct = initialize_spring_plot(num_zigs, w);
    box = rectangle('Position', []); % [x, y, w, h]

end

%% Will generate plots of springs

function spring_plotting_example()
    
    num_zigs = 5;
    w = .1;
    
    hold on;
    spring_plot_struct = initialize_spring_plot(num_zigs,w);
    axis equal; axis square;
    axis([-3, 3, -3, 3]);
    
    for theta = linspace(0, 6*pi, 1000)
        P1 = [0.5; 0.5];
        P2 = 2*[cos(theta); sin(theta)];
        update_spring_plot(spring_plot_struct, P1, P2)
        drawnow;
    end
end

%% Updates spring plotting object so that spring is plotted with ends 
% located at points P1 and P2

function update_spring_plot(spring_plot_struct, P1, P2)
    
    dP = P2 - P1;
    R = [dP(1), -dP(2)/norm(dP); dP(2), dP(1)/norm(dP)];
    
    plot_pts = R*spring_plot_struct.zig_zag;
    
    set(spring_plot_struct.line_plot,...
        'xdata',plot_pts(1,:)+P1(1),...
        'ydata',plot_pts(2,:)+P1(2));
    set(spring_plot_struct.point_plot,...
        'xdata',[P1(1),P2(1)],...
        'ydata',[P1(2),P2(2)]);

end

%% Create a struct containing plotting info for a single spring

% INPUTS:
%   num_zigs: number of zig zags in spring drawing
%   w: width of the spring drawing

function spring_plot_struct = initialize_spring_plot(num_zigs, w)
    
    spring_plot_struct = struct();
    
    zig_ending = [0.25, 0.75, 1; -1, 1, 0];
    zig_zag = zeros(2, 3 + 3*num_zigs);
    zig_zag(:, 1) = [-0.5; 0];
    zig_zag(:, end) = [num_zigs+0.5; 0];

    for n = 0:(num_zigs-1)
        zig_zag(:, (3+3*n):2+3*(n+1)) = zig_ending + [n, n, n; 0, 0, 0];
    end

    zig_zag(1,:) = (zig_zag(1,:)-zig_zag(1,1)) / (zig_zag(1,end)-zig_zag(1,1));
    zig_zag(2,:)= zig_zag(2,:)*w;
    
    spring_plot_struct.zig_zag = zig_zag;
    spring_plot_struct.line_plot = plot(0, 0, 'k', 'linewidth', 2);
    spring_plot_struct.point_plot = plot(0, 0,'ro', 'markerfacecolor', 'r', 'markersize', 7);

end

%% Initialize Orion Param System

function box_params = get_Orion_params()

    LW = 10; LH = 1; LG = 3;
    m = 1; Ic = (1/12)*(LH^2 + LW^2);
    g = 1; k = 20; k_list = [0.5*k, 0.5*k, 2*k, 5*k];
    l0 = 1.5*LG;
    
    Pbl_box = [-LW; -LH]/2;
    Pbr_box = [LW; -LH]/2;
    Ptl_box = [-LW; LH]/2;
    Ptr_box = [LW; LH]/2;
    
    boundary_pts = [Pbl_box, Pbr_box, Ptr_box, Ptl_box, Pbl_box];
    Pbl1_world = Pbl_box + [-LG; -LG];
    Pbl2_world = Pbl_box + [LG; -LG];
    
    Pbr1_world = Pbr_box + [0; -l0];
    Pbr2_world = Pbr_box + [l0; 0];
    P_world = [Pbl1_world, Pbl2_world, Pbr1_world, Pbr2_world];
    P_box = [Pbl_box, Pbl_box, Pbr_box, Pbr_box];
    
    % define system parameters
    box_params = struct();
    box_params.m = m;
    box_params.I = Ic;
    box_params.g = g;
    box_params.k_list = k_list;
    box_params.l0_list = l0*ones(size(P_world,2));
    box_params.P_world = P_world;
    box_params.P_box = P_box;
    box_params.boundary_pts = boundary_pts;

end