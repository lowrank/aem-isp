function [ val ] = neumann( x )
% A Neumann boundary function

val = (sin(10 * x(1,:)) .* (x(1,:) == 1)) + (0 * (x(1,:) == 0)) + ...
    (0 * (x(2,:) == 0)) + (0* (x(2,:) == 1));

end

