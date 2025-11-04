%% Computes the force exerted by the spring at one of its ends

% INPUTS:
%   PA: the position of the first end of the spring
%   PB: the position of the second end of the spring
%   k: the stiffness of the spring
%   l0: the natural length of the spring
% OUTPUTS:
%   F: the force exerted by the spring at end B

function F = compute_spring_force(k, l0, PA, PB)

% current length of the spring
l = norm(PA - PB);

% unit vector pointing from PA to PB
e_s = (PA - PB)/l;

% Force exerted by spring at point B
F = -k*(l - l0)*e_s;

end