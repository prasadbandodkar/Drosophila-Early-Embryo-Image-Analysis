function [H, I, bcoor] = find_peaks_guided(t, nPeaks, s0, cfg)
% FIND_PEAKS_GUIDED  Self-contained peak finder guided by canonical locations.
%
%   [H, I, bcoor] = find_peaks_guided(t, nPeaks, s0)
%
%   Inputs:
%     t       - 1D signal (row or column vector)
%     nPeaks  - number of peaks to return
%     s0      - canonical peak locations in normalised [0,1] coordinates,
%               length nPeaks
%
%   Outputs:
%     H     - [1 x nPeaks] peak heights  (NaN where not found)
%     I     - [1 x nPeaks] peak indices  (NaN where not found)
%     bcoor - [nPeaks x 2] boundary indices [left, right] for each peak
%
%   No toolbox dependencies — local maxima and prominences are computed
%   internally. All tuning parameters are derived automatically from the
%   signal and s0.

arguments
    t      (1,:) double
    nPeaks (1,1) double {mustBePositive, mustBeInteger}
    s0     (1,:) double
    cfg    struct
end

t = t(:)';                          % enforce row
n = length(t);
s = linspace(0, 1, n);

% ---- auto parameters -------------------------------------------------
sigRange         = max(t) - min(t);
minPeakProm      = 0.05 * sigRange;
if nPeaks > 1
    tolerance    = 0.5 * min(diff(s0));
else
    tolerance    = 0.25;
end

% ---- step 1: find all local maxima ----------------------------------
% interior: strictly greater than both neighbors
% endpoints: strictly greater than their single neighbor
isMax = [t(1) > t(2), ...
         t(2:end-1) > t(1:end-2) & t(2:end-1) > t(3:end), ...
         t(end) > t(end-1)];
Icand = find(isMax);
Hcand = t(Icand);

% ---- step 2: compute prominence for each candidate ------------------
% prominence = height - max(left_base, right_base)
% base = lowest trough between the peak and any higher peak
% endpoints have no data on one side, so that base is -Inf
prom = zeros(size(Icand));
for ci = 1:length(Icand)
    pk = Icand(ci);
    h  = Hcand(ci);

    % left base: lowest value between pk and the nearest higher peak to left
    left_higher = Icand(Icand < pk & Hcand > h);
    if pk == 1
        left_base = -Inf;                          % no data to the left
    elseif isempty(left_higher)
        left_base = min(t(1:pk));
    else
        left_base = min(t(left_higher(end):pk));
    end

    % right base: lowest value between pk and nearest higher peak to right
    right_higher = Icand(Icand > pk & Hcand > h);
    if pk == n
        right_base = -Inf;                         % no data to the right
    elseif isempty(right_higher)
        right_base = min(t(pk:end));
    else
        right_base = min(t(pk:right_higher(1)));
    end

    prom(ci) = h - max(left_base, right_base);
end

% ---- step 3: filter by MinPeakProminence ----------------------------
keep  = prom >= minPeakProm;
Icand = Icand(keep);
Hcand = Hcand(keep);

% ---- step 4: assign best candidate per s0(k) ------------------------
H = nan(1, nPeaks);
I = nan(1, nPeaks);

for k = 1:nPeaks
    d      = abs(s(Icand) - s0(k));
    within = d <= tolerance;

    % Windowed maximum: always the global max within the tolerance window,
    % regardless of local-max / prominence filtering.  This is the prior-
    % guided fallback and ensures endpoints or monotone ramps are handled.
    lo  = max(1, round((s0(k) - tolerance) * n));
    hi  = min(n, round((s0(k) + tolerance) * n));
    seg = t(lo:hi);
    if isempty(seg)
        [~, wI] = min(abs(s - s0(k)));
    else
        [~, rel] = max(seg);
        wI = lo + rel - 1;
    end
    wH = t(wI);

    if any(within)
        % Merge filtered candidates with the windowed max, pick the tallest.
        allH = [Hcand(within), wH];
        allI = [Icand(within), wI];
        [H(k), best] = max(allH);
        I(k) = allI(best);
    else
        H(k) = wH;
        I(k) = wI;
    end
end

% ---- step 5: compute peak boundaries (bcoor) ------------------------
% bcoor(k,:) = [left_idx, right_idx] of the fitting window for peak k.
% Uses 10% of peak height as threshold for edge detection.
bh    = 0.10;
bcoor = zeros(nPeaks, 2);

% Left boundary of first peak: last crossing of bh*H(1) before I(1)
cross = find((t - bh*H(1)) .* (t([2:end, 1]) - bh*H(1)) < 0);
jf    = find(cross < I(1));
if isempty(jf)
    bcoor(1, 1) = 1;
else
    jf = jf(end);
    loc = 0.5 * (s(cross(jf)) + s(cross(jf) + 1));
    [~, bcoor(1, 1)] = min(abs(s - loc));
end

% Boundaries between adjacent peaks: index of minimum between consecutive peaks
if nPeaks > 1
    for k = 2:nPeaks
        t_seg = t(I(k-1):I(k));
        [~, i0] = min(t_seg);
        bcoor(k-1, 2) = I(k-1) + i0 - 1;
        bcoor(k,   1) = bcoor(k-1, 2);
    end
end

% Right boundary of last peak: first crossing of bh*H(end) after I(end)
cross = find((t - bh*H(end)) .* (t([2:end, 1]) - bh*H(end)) < 0);
jf    = find(cross > I(end));
if isempty(jf) || isempty(cross)
    bcoor(end, 2) = n;
else
    jf = jf(1);
    if cross(jf) + 1 <= n
        loc = 0.5 * (s(cross(jf)) + s(cross(jf) + 1));
        [~, bcoor(end, 2)] = min(abs(s - loc));
    else
        bcoor(end, 2) = n;
    end
end


% Finally, plot the peaks if cfg.debug mode is true
if cfg.debug
    figure('Visible', 'off');
    plot(s, t, 'k-');
    hold on;
    plot(s(I), H, 'ro', 'MarkerFaceColor', 'r');
    plot(s(bcoor(:,1)), t(bcoor(:,1)), 'go', 'MarkerFaceColor', 'g');
    plot(s(bcoor(:,2)), t(bcoor(:,2)), 'bo', 'MarkerFaceColor', 'b');
    [~, Is0] = arrayfun(@(x) min(abs(s - x)), s0);
    plot(s0, t(Is0), 'm^', 'MarkerFaceColor', 'm');
    hold off;
    if strcmp(cfg.image_type, 'cross-section')
        legend('Signal', 'Peaks', 'Ventral Border', 'Dorsal Border', 'Canonical locations', 'Location', 'best');
    else
        legend('Signal', 'Peaks', 'Anterior Border', 'Posterior Border', 'Canonical locations', 'Location', 'best');
    end
    title(cfg.plot_title);
    xlabel('Position');
    ylabel('Intensity');
    saveas(gcf, fullfile(cfg.path_data, ['find_peaks_', cfg.plot_filename]));
end

end
