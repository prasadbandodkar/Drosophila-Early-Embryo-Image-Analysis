function [data, meta] = fit_peaks(t, channelname, nPeaks, H, I, bcoor, cfg)


hh = 0.5;
ns = length(t);
s  = linspace(0, 1, ns)';
t  = t(:);          % ensure column vector

tFit = 0;
for i = 1:nPeaks

    % Extract windowed peak signal
    ti = zeros(ns, 1);
    ti(bcoor(i,1):bcoor(i,2)) = t(bcoor(i,1):bcoor(i,2));

    % Build fittype for this peak
    if nPeaks > 1
        genename = [channelname, '_', num2str(i)];
    else
        genename = channelname;
    end
    ii         = num2str(i);
    params     = ['A', ii, 'a', ',delt', ii, 'a', ',x0', ii, 'a'];
    fittypestr  = ['genefit(''', genename, ''',x,B,', params, ')'];
    f           = fittype(fittypestr);
    [cfun, cvals, cint, cint68, rsquare] = fitelephant(ti, s, H(i), I(i), f);

    % Canonical gene borders at half-height
    [xA, xP, x_offset] = canonicalgeneborders(genename, hh);
    if isnan(xA), xA = 0; end
    if isnan(xP), xP = 1; end

    % Unpack fit coefficients
    A(i)  = cvals(1);
    B(i)  = cvals(2);
    delt  = cvals(3);
    % For ventral/dorsal genes (x_offset == 0 or 1), genefit fixes x0 at
    % x_offset and ignores the x0 coefficient — mirror that here.
    if x_offset == 0 || x_offset == 1
        x0 = x_offset;
    else
        x0 = cvals(4);
    end
    dcint = 0.5 * diff(cint68);
    ddelt = dcint(3);
    dx0   = dcint(4);

    % Derived spatial quantities
    data.sPeaks(i)     = x0;
    data.LeftBorder(i)  = (xA - x_offset) * delt + x0;
    data.RightBorder(i) = (xP - x_offset) * delt + x0;
    data.Width(i)       = data.RightBorder(i) - data.LeftBorder(i);
    data.dsPeaks(i)     = dx0;
    data.dLeftBorder(i) = sqrt((xA - x_offset)^2 * ddelt^2 + dx0^2);
    data.dRightBorder(i)= sqrt((xP - x_offset)^2 * ddelt^2 + dx0^2);
    data.dWidth(i)      = data.dRightBorder(i) - data.dLeftBorder(i);

    % Fit metadata
    meta.Delt(i)    = delt;
    meta.Rsquare(i) = rsquare;
    meta.Cint{i}    = cint;
    meta.Cint68{i}  = cint68;

    % Accumulate fit profile
    tFit = tFit + cfun(s);
end

meta.tFit = tFit;

% Score: RMS residual over the region spanning first LeftBorder to last RightBorder
[~, iLeft]  = min(abs(s - data.LeftBorder(1)));
[~, iRight] = min(abs(s - data.RightBorder(nPeaks)));
i0          = iLeft:iRight;
score       = sqrt(sum((tFit(i0) - t(i0)).^2)) / length(t(i0));
data.score  = score;


if cfg.debug
    figure('Visible', 'off');
    plot(s, t, 'k-');
    hold on;
    plot(s, tFit, 'r-', 'LineWidth', 1.5);
    plot(s(I), H, 'ro', 'MarkerFaceColor', 'r');
    if strcmp(cfg.image_type, 'cross-section')
        leftBorderName = 'Ventral Border';
        rightBorderName = 'Dorsal Border';
    else
        leftBorderName = 'Anterior Border';
        rightBorderName = 'Posterior Border';
    end
    for i = 1:nPeaks
        xline(data.LeftBorder(i), 'g-', sprintf('%s%d', leftBorderName, i));
        xline(data.RightBorder(i), 'b-', sprintf('%s%d', rightBorderName, i));
    end
    hold off;
    legend('Signal', 'Fit', 'Peaks', leftBorderName, rightBorderName, 'Location', 'best');
    title(cfg.plot_title);
    xlabel('Position');
    ylabel('Intensity');
    saveas(gcf, fullfile(cfg.path_data, ['fit_peaks_', cfg.plot_filename]));
end

end
