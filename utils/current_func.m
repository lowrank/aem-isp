function [ ret ] = current_func( x )


% ret1 = 0.1 + 0.02 * cos(2 * x(1,:)) .* cos(2 * x(2,:));
% ret2 = 0.1 + 0.02 * cos(4 * (x(1,:).^2 + x(2,:).^2));
% output is 2 x N matrix
ret = [derenzo(flipud(1.25 * (x-0.5) + 0.5), 0.0, 16); ...
    derenzo(1.45 * (x-0.5) + 0.5, 0.0, 8)]; % vanishes on boundary



end

