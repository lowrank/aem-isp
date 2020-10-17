function [ ret ] = shepp_logan( x, N, scale, shift)

%ABSORPTIONFF Absorption coef for fluorescence material.
%   x coordinate, vectorized.
%   ret ~0.2
if nargin < 2
    N = 500;
    scale = 0.1;
    shift = 0.05;
end

[X,Y] =meshgrid(linspace(0,1, N));
X = X(:);
Y = Y(:);


Bone = 6.25e-5;
Grey = 3.334e-3;
CSF  = 15.38e-3;
White = 1.428e-3;
Background = 4.3478e-3;
%  Background = 0.1;


coeff = [  (Bone - Background)   .69   .92    0     0     0   
        (Grey - Bone)  .6624 .8740   0  -.0184   0
        (CSF - Grey)  .1100 .3100  .22    0    -18
        (CSF - Grey)  .1600 .4100 -.22    0     18
        (White - Grey)  .2100 .2500   0    .35    0
         (CSF - Grey)  .0460 .0460   0    .1     0
         (CSF - Grey)  .0460 .0460   0   -.1     0
         (CSF - Grey)  .0460 .0230 -.08  -.605   0 
        (CSF - Grey)  .0230 .0230   0   -.606   0
         (CSF - Grey)  .0230 .0460  .06  -.605   0   ];


ph = phantom(coeff, N);

% ph = imgaussian(ph, 20);

F = scatteredInterpolant(X,Y,ph(:));
ret = F(x(1,:), x(2,:)) * scale + shift;


end

