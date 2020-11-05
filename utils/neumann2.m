function [ val ] = neumann2( x )
% Another Neumann boundary function

theta = pi ;
val = (cos(theta)* (x(1,:) == 0.9)) + (-cos(theta) .* (x(1,:) == 0.1)) + ...
    (-sin(theta) * (x(2,:) == 0)) + (sin(theta) * (x(2,:) == 1));

end
