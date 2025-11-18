%Implementation of finite difference approximation
%for Jacobian of multidimensional function
%INPUTS:
%fun: the mathetmatical function we want to differentiate
%x: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%J: approximation of Jacobian of fun at x
function J = approximate_jacobian(fun, x)
    m = length(fun(x));
    J = zeros(m, length(x));
    e = [];
    h = 1e-6;
    for i = 1:length(x)
        e = zeros(length(x), 1);
        e(i) = 1;
        
        x_plus = x + h*e;
        x_minus = x - h*e;

        new_column = (fun(x_plus) - fun(x_minus)) / (2*h);
        J(:,i) = new_column;
    end
end