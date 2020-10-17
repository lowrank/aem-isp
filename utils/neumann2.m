function [ val ] = neumann2( x )
% Another Neumann boundary function

theta = 1.5 * pi/4;
val = (cos(theta)* (x(1,:) == 1)) + (-cos(theta) .* (x(1,:) == 0)) + ...
    (-sin(theta) * (x(2,:) == 0)) + (sin(theta) * (x(2,:) == 1));

end
