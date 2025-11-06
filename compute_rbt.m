%% Computes the rigid body transformation that maps a set of points in the 
% box-frame to their corresponding world-frame coordinates

% INPUTS:
%   x: the x position of the centroid of the box
%   y: the y position of the centroid of the box
%   theta: the orientation of the box
%   Plist_box: a 2 x n matrix of points in the box frame
% OUTPUTS:
%   Plist_world: a 2 x n matrix of points describing the world-frame 
%           coordinates of the points in Plist_box [x_vals; y_vals]

function Plist_world = compute_rbt(x, y, theta, Plist_box)
    
    rc_world = [x; y]; % position of centroid as vector

    % construct CCW rotation matrix with theta
    ccw_rot_mat = [cos(theta), -sin(theta); sin(theta), cos(theta);]; 
    
    % Rotate: qi_world = [CCW rotation matrix]*[qi_body]
    qi_world = ccw_rot_mat*Plist_box;
    % Translate: ri_world = rc_world + qi_world
    Plist_world = rc_world + qi_world;

end