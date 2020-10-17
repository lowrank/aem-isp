function [ val ] = neumann( x )
% A Neumann boundary function

theta = 1* pi/4;
% cos theta, sin theta

val = (cos(theta)* (x(1,:) == 1)) + (-cos(theta) .* (x(1,:) == 0)) + ...
    (-sin(theta) * (x(2,:) == 0)) + (sin(theta) * (x(2,:) == 1));

end
