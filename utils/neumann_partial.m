function [ val ] = neumann_partial( x, mask )
% A Neumann boundary function

theta =0;
% cos theta, sin theta

val = (cos(theta)* (x(1,:) == 0.9)) + (-cos(theta) .* (x(1,:) == 0.1)) + ...
    (-sin(theta) * (x(2,:) == 0)) + (sin(theta) * (x(2,:) == 1));

val = val .* mask(x);

end
