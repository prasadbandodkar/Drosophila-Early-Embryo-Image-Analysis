function [xp_in, yp_in, xp_in2, yp_in2] = find_normals(xp, yp, options)
% FIND_NORMALS  Find points locally normal to a boundary curve.
%
%   [xp_in, yp_in, xp_in2, yp_in2] = find_normals(xp, yp)
%   [xp_in, yp_in, xp_in2, yp_in2] = find_normals(xp, yp, Yhatmax=60, dthetamax=45)
%
%   Given boundary points (xp, yp), returns two offset contours shifted
%   inward by Yhatmax pixels along the local normal direction.

arguments
    xp         (:,1) double
    yp         (:,1) double
    options.Yhatmax   (1,1) double = 60
    options.dthetamax (1,1) double = 45
end

Yhatmax = options.Yhatmax;

ns = length(xp) - 1;

% Tangent vectors (centered difference, wrapping around endpoints)
a1 = xp(2:end) - [xp(end-1); xp(1:end-2)];
a2 = yp(2:end) - [yp(end-1); yp(1:end-2)];

% Unit normals: rotate tangent 90° and normalize
% Normal to (a1,a2) is (-a2, a1); normalize by magnitude
mag = sqrt(a1.^2 + a2.^2);
b1  = -a2 ./ mag;
b2  =  a1 ./ mag;

% Offset inner contour 1: aligned with xp(1:end-1), closed at end
xp_in  = [xp(1:end-1) + Yhatmax * b1; 0];
yp_in  = [yp(1:end-1) + Yhatmax * b2; 0];
xp_in(end) = xp_in(1);
yp_in(end) = yp_in(1);

% Offset inner contour 2: aligned with xp(2:end), closed at start
xp_in2    = zeros(ns+1, 1);
yp_in2    = zeros(ns+1, 1);
xp_in2(2:end) = xp(2:end) + Yhatmax * b1;
yp_in2(2:end) = yp(2:end) + Yhatmax * b2;
xp_in2(1) = xp_in(end);
yp_in2(1) = yp_in(end);

end
