function [anterior, posterior] = find_poles(xp, yp, phi)
% find_poles - Find the anterior and posterior poles of an embryo periphery.
%
% Fits the periphery to an ellipse, rotates so the AP axis is horizontal,
% finds the two extreme endpoints, and computes curvature to identify poles.
%
% Inputs:
%   xp, yp   - x and y coordinates of the embryo periphery
%
% Optional name-value inputs:
%   phi      - ellipse orientation angle (radians); if omitted, computed
%              automatically via ellipse_fit
%
% Outputs:
%   anterior  - structs.pole for the anterior pole (smaller radius of curvature)
%   posterior - structs.pole for the posterior pole

arguments
    xp
    yp
    phi (1,1) double = NaN
end

% Fit ellipse (or use provided phi) and rotate so the major axis is horizontal
if isnan(phi)
    phi = functions.ellipse_fit(xp, yp);
else
    phi = phi;
end
if phi > pi/2
    phi = phi - pi;
end
[xp2, yp2] = rotate_points(pi - phi, xp, yp);

% Pad point arrays to handle boundary wrapping
n = length(xp2);
x = [xp2; xp2(2:end-1); xp2];
y = [yp2; yp2(2:end-1); yp2];

% Find the two pole endpoints (min/max x in rotated frame)
nCurve = round(n / 25);
[~, i1] = min(xp2);  i1 = i1 + n;
[~, i2] = max(xp2);  i2 = i2 + n;

xPole1 = x(i1-nCurve:i1+nCurve);  yPole1 = y(i1-nCurve:i1+nCurve);
xPole2 = x(i2-nCurve:i2+nCurve);  yPole2 = y(i2-nCurve:i2+nCurve);

% Compute curvature at each pole
[xP1, rad1] = find_curvature(xPole1, yPole1);
[xP2, rad2] = find_curvature(xPole2, yPole2);

% Match curvature x-positions back to original indices
[~, i1] = min(abs(xP1 - xp2));
[~, i2] = min(abs(xP2 - xp2));

% Arclength coordinate at each pole index
ds    = sqrt(diff(xp).^2 + diff(yp).^2);
ss    = [0; cumsum(ds)] / sum(ds) * 2 - 1;  % [-1, 1]

% Anterior pole has smaller radius of curvature (sharper tip)
if rad1 < rad2
    anterior  = structs.pole(idx=i1, s=ss(i1), xp=xp(i1), yp=yp(i1), rad=rad1);
    posterior = structs.pole(idx=i2, s=ss(i2), xp=xp(i2), yp=yp(i2), rad=rad2);
else
    anterior  = structs.pole(idx=i2, s=ss(i2), xp=xp(i2), yp=yp(i2), rad=rad2);
    posterior = structs.pole(idx=i1, s=ss(i1), xp=xp(i1), yp=yp(i1), rad=rad1);
end

end


% --- Sub-function: radius of curvature at the tip of a pole segment ------
function [xend, Rmin] = find_curvature(xpend, ypend)
    y       = linspace(min(ypend), max(ypend), 100);
    p       = polyfit(ypend, xpend, 2);
    x       = polyval(p, y);
    radCurv = abs((1 + (2*p(1)*y + p(2)).^2).^1.5 / (2*p(1)));
    [Rmin, iend] = min(radCurv);
    xend = x(iend);
end


% --- Sub-function: rotate 2D points by angle theta -----------------------
function [x_out, y_out] = rotate_points(theta, x, y)
    x_out = x*cos(theta) - y*sin(theta);
    y_out = y*cos(theta) + x*sin(theta);
end
