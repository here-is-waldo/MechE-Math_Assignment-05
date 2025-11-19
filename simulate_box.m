
function simulate_box()
    
    % get params
    box_params = get_Orion_params();
    LW = 10; LH = 1; % from orion params, may want to add to box_params
    
    % load the system parameters into the rate function
    % via an anonymous function
    my_rate_func = @(t_in, V_in) box_rate_func(t_in, V_in, box_params);
    
    % more initial conditions
    x0 = 0; y0 = 0; theta0 = 0;
    vx0 = 5; vy0 = 5; vtheta0 = pi/6;
    V0 = [x0; y0; theta0; vx0; vy0; vtheta0];
    tspan = [0, 3];

    V_eq = find_equilibrium(box_params, V0); % Compute equilibrium state
    [A, J_approx] = linearize_system(box_params, V_eq);
    Q = -A(4:6, 1:3);
    [U_mode, omega_n] = modal_analysis(Q, V_eq, box_params);
    V0 = V_eq;
    clf;
%% 

    % run the integration (using ode45 cuz my integrator isn't reliable)
    % THIS STEP IS BIG WEIRD 
    % [tlist, Vlist] = ode45(my_rate_func, tspan, V0);
    h_ref = 0.01; BT_struct = get_BT("Dormand Prince");
    [tlist, Vlist, ~, ~] = ARI_explicit_RK_fixed_step_integration(my_rate_func, tspan, V0, h_ref, BT_struct);


    %%%%% ANIMATION %%%%%

    % File set up:
    keep_vid = false; % record and store a video bool
    if keep_vid == true

        % define location and filename where video will be stored
        mypath1 = 'C:\Users\lodio\OneDrive - Olin College of Engineering\Desktop\';
        mypath2 = 'Classes\Junior Fall\Orion Math\Assignment-05\MechE-Math_Assignment-05';
        fname = 'debug2.avi';
        input_fname = [mypath1, mypath2, fname];

        % create a videowriter, which will write frames to the animation file
        writerObj = VideoWriter(input_fname);
        % must call open before writing any frames
        open(writerObj);
    end
    
    % Input set up:
    % initialize vars
    num_zigs = 8; w = -.1;
    num_springs = length(box_params.P_world);
    tdiff = [0; diff(tlist)];
    x = V0(1); y = V0(2); theta = V0(3);

    % initialize figure
    fig1 = figure('Name','Box Animation','NumberTitle','off');
    clf(fig1);                  % clear but do NOT steal figures created elsewhere
    ax = axes('Parent', fig1);  % create dedicated axes
    axis(ax, 'equal');
    hold(ax, 'on');

    % initialize plots 
    all_spring_plots = cell(num_springs, 1);
    all_spring_plots{1} = initialize_spring_plot(num_zigs, w); hold on;
    for k = 2:num_springs
        all_spring_plots{k} = initialize_spring_plot(num_zigs, w);
    end
    box_plot_struct = initialize_box_plot(x, y, theta, LH, LW, box_params);

    % Animation:
    % loop and plot each timestep
    for i = 1:length(tlist)
        disp(i)
        % delay according to timestep
        pause(tdiff(i));

        % index for current state
        x = Vlist(i,1); y = Vlist(i,2); theta = Vlist(i,3);

        % update plots
        update_box_plot(box_plot_struct, x, y, theta, box_params);
        % Note: update_box_plot() updates box_params.boundary_pts and
        % box_params.P_box
        
        % used updated P_box to plot springs
        P1_list = compute_rbt(x, y, theta, box_params.P_box);
        P2_list = box_params.P_world;
        for j = 1:length(P1_list)
            update_spring_plot(all_spring_plots{j}, P1_list(:,j), P2_list(:,j));
        end
        
        % redraw
        drawnow;
        axis equal;

        if keep_vid == true
            % capture a frame (what is currently plotted)
            current_frame = getframe(fig1);
            % write the frame to the video
            writeVideo(writerObj, current_frame);
        end
    end

    % Clean up: close writer object if recording
    if keep_vid == true
        close(writerObj)
    end

end

%% Updates spring plotting object so that spring is plotted with ends 
% located at points P1 and P2

function update_spring_plot(spring_plot_struct, P1, P2)
    
    dP = P2 - P1;
    R = [dP(1), -dP(2)/norm(dP); dP(2), dP(1)/norm(dP)];
    
    plot_pts = R*spring_plot_struct.zig_zag;

    set(spring_plot_struct.line_plot, ...
        'xdata', plot_pts(1,:)+P1(1), ...
        'ydata', plot_pts(2,:)+P1(2));
    set(spring_plot_struct.point_plot, ...
        'xdata', [P1(1), P2(1)], ...
        'ydata', [P1(2), P2(2)]);

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
    spring_plot_struct.line_plot = plot(0, 0, 'k', 'linewidth', 2); hold on;
    spring_plot_struct.point_plot = plot(0, 0,'ro', 'markerfacecolor', 'r', 'markersize', 7);

end

%% Updates box plotting object so that the box is plotted in accordance with
% state x, y, theta

function update_box_plot(box_plot_struct, x, y, theta, box_params)
    
    % calc & update boundary points in box frame
    % box_params.boundary_pts = [x-w, x+w, x+w, x-w, x-w; ...
    %                            y-h, y-h, y+h, y+h, y-h;]/2;
    % update P_box points (SYSTEM SPECIFIC)
    Pbl_box = box_params.boundary_pts(:,1);
    Pbr_box = box_params.boundary_pts(:,2);
    box_params.P_box = [Pbl_box, Pbl_box, Pbr_box, Pbr_box];
    
    % convert to world frame
    box_pts_world = compute_rbt(x, y, theta, box_params.boundary_pts);
    
    % update data via set()
    set(box_plot_struct.line_plot, ...
        'xdata', box_pts_world(1, :), ...
        'ydata', box_pts_world(2, :));
    set(box_plot_struct.point_plot, ...
        'xdata', box_pts_world(1, :), ...
        'ydata', box_pts_world(2, :));

end

%% Create a struct containing plotting info for the box given (x, y, z)

% INPUTS:
%   x: current x position of the box
%   y: current y position of the box
%   theta: current orientation of the box
%   h: height of the box (oriented in the direction of theta)
%   w: width of the box (oriented perpendicular to theta)
% OUTPUTS
%   box_plot_struct: a struct containing the line plot and the point
%                   plot objects for the box

function box_plot_struct = initialize_box_plot(x, y, theta, h, w, box_params)

    % calc & update boundary points in box frame
    box_params.boundary_pts = [x-w, x+w, x+w, x-w, x-w; ...
                               y-h, y-h, y+h, y+h, y-h;]/2;
    % update P_box points (SYSTEM SPECIFIC)
    Pbl_box = box_params.boundary_pts(:,1);
    Pbr_box = box_params.boundary_pts(:,2);
    box_params.P_box = [Pbl_box, Pbl_box, Pbr_box, Pbr_box];

    % convert to world frame
    box_pts_world = compute_rbt(x, y, theta, box_params.boundary_pts);

    % store in vars for convenience
    x_pts = box_pts_world(1,:);
    y_pts = box_pts_world(2,:);

    % plot lines and points
    box_plot_struct = struct();
    box_plot_struct.line_plot = plot(x_pts, y_pts, 'k-', 'linewidth', 2); hold on;
    axis equal;
    box_plot_struct.point_plot = plot(x_pts, y_pts, 'ro', 'markerfacecolor', 'r', 'markersize', 7);

end

%% Initialize Orion Param System

function box_params = get_Orion_params()

    LW = 10; LH = 1; LG = 3;
    m = 1; Ic = (1/12)*(LH^2 + LW^2);
    g = 1; k = 5; k_list = [0.5*k, 0.5*k, 2*k, 5*k];
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
    
    % assign system parameters
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