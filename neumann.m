function [ val ] = neumann( x )
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

val = (sin(10 * x(1,:)) .* (x(1,:) == 1)) + (0 * (x(1,:) == 0)) + ...
    (0 * (x(2,:) == 0)) + (0* (x(2,:) == 1));

end

