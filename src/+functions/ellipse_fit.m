function phi = ellipse_fit(x, y)
% Fit ellipse to points (x,y) and return orientation angle phi (radians).

x = x(:);
y = y(:);

% Least-squares fit: solve for conic coefficients (a=1 fixed)
M = [2*x.*y, y.^2, 2*x, 2*y, ones(size(x))];
e = M \ (-x.^2);

% Conic coefficients
a = 1;
b = e(1);
c = e(2);

% Orientation angle (Wolfram Mathworld formula)
if a < c
    phi = 0.5 * acot((a - c) / (2*b));
else
    phi = pi/2 + 0.5 * acot((a - c) / (2*b));
end
