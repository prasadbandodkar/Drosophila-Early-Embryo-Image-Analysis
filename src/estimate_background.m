function [bg,bgsig] = estimate_background(I, dtype)
% Estimate background intensity level for a single image channel.
%
% We assume the background intensity level (i.e., the true black level) is
% equal to the mode of pixel intensities in the image, estimated by finding
% the peak of the histogram and fitting a local Gaussian to it.
%
% Inputs:
%   I     - 2D image array (single channel)
%   dtype - pixel type string, e.g. 'uint8', 'uint16'
%
% Outputs:
%   bg    - estimated background intensity (mode of histogram)
%   bgsig - estimated background noise (std dev of Gaussian fit)

switch dtype
    case 'uint8'
        nbins = 256;
    case 'uint16'
        nbins = 65536;
    otherwise
        error('estimate_background: unsupported dtype "%s"', dtype);
end

% Get histogram of pixel intensities
edges      = 0:nbins;                        % nbins+1 edges → nbins bins
[n, edges] = histcounts(double(I(:)), edges);
x          = edges(1:end-1) + 0.5;          % bin centres

% Find the peak of the histogram
[A, k] = max(n(1:end-1));                   % ignore last (saturated) bin
bg      = x(k);

% Fit a Gaussian to the peak
X     = x(max(k-4,1) : k+4);
Y     = n(max(k-4,1) : k+4);
SIG   = (X - bg) ./ sqrt(-2*log(Y/A));
bgsig = sqrt(mean(SIG'.^2, 'omitnan'));

end
