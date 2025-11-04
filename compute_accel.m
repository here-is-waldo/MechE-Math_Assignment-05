%% Computes the linear and angular acceleration of the box given its 
% current position and orientation

% INPUTS:
%   x: current x position of the box
%   y: current y position of the box
%   theta: current orientation of the box
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
% OUTPUTS
%   ax: x acceleration of the box
%   ay: y acceleration of the box
%   atheta: angular acceleration of the box

function [ax, ay, atheta] = compute_accel(x, y, theta, box_params)

    % EQUATIONS:
    % theta'' = 1/I_c * sum(qi_world x Fi_world)
    % m*r_centroid = -mg j^ + sum(Fi_world) (the spring forces)

    % unpack vars
    m = box_params.m;
    I_c = box_params.I;
    g = box_params.g;
    r_c = [x; y]; % pos vector of centroid

    % intermediate computations
    sum_Fi = 0; % sum of all spring forces
    sum_qi_x_Fi = 0; % sum of all qi x Fi

    % iterate through each spring
    for i = 1:size(box_params.P_box, 2)

        % unpack vars
        k = box_params.k_list(i);
        l0 = box_params.l0_list(i);
        
        % unpack and transform spring mounting positions
        PB_world = box_params.P_world(:, i);
        PA_box = box_params.P_box(:, i); % in box reference frame
        PA_world = compute_rbt(x, y, theta, PA_box); % now world frame

        % find force and distance from centroid
        Fi = compute_spring_force(k, l0, PA_world, PB_world);
        qi = r_c - PA_world;

        % add to total
        sum_Fi = sum_Fi + Fi;
        sum_qi_x_Fi = sum_qi_x_Fi + cross(qi, Fi);
    end

    % Linear acceleration
    a_vec = [0, -g] + (1/m)*sum_Fi; 
    ax = a_vec(1); ay = a_vec(2);

    % Angular acceleration
    atheta = (1/I_c)*sum_qi_x_Fi; 

end