function [ ret ] = current_func( x )




% ret1 = 0.1 + 0.02 * cos(2 * x(1,:)) .* cos(2 * x(2,:));
% ret2 = 0.1 + 0.02 * cos(4 * (x(1,:).^2 + x(2,:).^2));
% output is 2 x N matrix
ret = [sin(2 * pi *  x(1,:)) .* sin(2 * pi *  x(2,:) ); ...
    sin(4 * pi *  x(1,:)) .* sin(4 * pi *  x(2,:) )]; % vanishes on boundary



end

