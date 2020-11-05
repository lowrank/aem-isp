function [ val ] = neumann( x )
% A Neumann boundary function

theta = 5 * pi/6;
% cos theta, sin theta

val = (cos(theta)* (x(1,:) == 0.9)) + (-cos(theta) .* (x(1,:) == 0.1)) + ...
    (-sin(theta) * (x(2,:) == 0)) + (sin(theta) * (x(2,:) == 1));

end
