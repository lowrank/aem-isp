function [ val ] = neumann2( x )
% Another Neumann boundary function

val = (0.2 * (x(1,:) == 1)) + (0.3 * (x(1,:) == 0)) + ...
    (1 * (x(2,:) == 0)) + (-0.5 * (x(2,:) == 1));

end
