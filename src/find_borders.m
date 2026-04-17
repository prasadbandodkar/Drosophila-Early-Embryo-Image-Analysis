function [xp, yp, I, seglabel] = find_borders(I, nt)
% find_borders  Find uniformly distributed coordinates of the embryo periphery.
%
%   [xp, yp, I, seglabel] = find_borders(I, nt)
%
%   Inputs:
%       I  - Image (uint16)
%       nt - Desired number of border points (forced even)
%
%   Outputs:
%       xp       - x-coordinates of the border (length nt+1, closed)
%       yp       - y-coordinates of the border (length nt+1, closed)
%       I        - Segmented image
%       seglabel - Label matrix from segmentation

% --- Segment and extract raw boundary ---
[I, seglabel] = segment_embryo_image(I);
boundary = bwboundaries(seglabel, 'noholes');
points   = boundary{1};
xp       = points(:, 2);
yp       = points(:, 1);

% --- Smooth periodically (wrap, movmean, unwrap) ---
N       = length(xp);
nsmooth = round(0.05 * N);
xp = smooth_periodic(xp, nsmooth);
yp = smooth_periodic(yp, nsmooth);

% --- Force nt even ---
if mod(nt, 2) ~= 0
    nt = nt + 1;
end

% --- Resample to nt uniformly arc-length spaced points (two passes) ---
for pass = 1:2
    ds   = sqrt(diff(xp).^2 + diff(yp).^2);
    ds(ds < eps) = eps;          % guard against zero-length segments
    s    = [0; cumsum(ds)];
    s_u  = linspace(0, s(end), nt + 1)';
    xp   = interp1(s, xp, s_u);
    yp   = interp1(s, yp, s_u);
    xp(end) = xp(1);             % close the curve exactly
    yp(end) = yp(1);
end

end

% -------------------------------------------------------------------------
function ysmooth = smooth_periodic(y, nsmooth)
% Smooth a periodic vector using a moving average.
    wasrow = isrow(y);
    y = y(:);
    y_wrapped = [y(end-nsmooth+1:end); y; y(1:nsmooth)];
    ysmooth   = movmean(y_wrapped, nsmooth);
    ysmooth   = ysmooth(nsmooth+1:end-nsmooth);
    if wasrow, ysmooth = ysmooth'; end
end