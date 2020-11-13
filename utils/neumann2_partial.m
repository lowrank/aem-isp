function [ val ] = neumann2_partial( x, mask )
% Another Neumann boundary function

theta = pi/2;
val = (cos(theta)* (x(1,:) == 0.9)) + (-cos(theta) .* (x(1,:) == 0.1)) + ...
   (-sin(theta) * (x(2,:) == 0)) + (sin(theta) * (x(2,:) == 1));

val = val .* mask(x);
end
