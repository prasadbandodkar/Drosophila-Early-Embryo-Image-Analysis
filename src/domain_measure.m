function [t,raw,s] = domain_measure(I,xp,yp,varargin)
%Finds domains of gene expression in a cross-sectioned embryo.
%
%function [t,raw,s] = domain_measure(I,xp,yp,varargin)
%

% This function takes a truecolor cross-sectional image of an embryo and
% finds the fluorescence intensity (gene expression) as a function of
% fractions of total circumference.
%
% Inputs:
%
% "I": can be either a uint8/uint16 grayscale image, or string of filename
% "xp,yp": the points along the perimeter of the embryo
% 
% Optional argument varargin can consist of the following things:
%	* "depthInsideImage":  the distance, in pixels into the embryo for
%       analysis
%		Default, 30 pixels.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%	* "ns": choice for number of bins in "s" (pseudoarclength) when
%		measuring the intensity around the periphery.  Must be an even
%		number.  Default, 300.
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% Outputs:
%
% "t": npts-by-nChannels array of smoothened fluorescent intensities.
% "raw": same as "t", but the raw (non-smoothened) data.
% "s": arclength, going from -1 to +1, with "npts" points.
% 
% PB Edit - 01/11/2020
% This function now finds the inner border, makes the border uniform (just 
% like the external border) and then finds the domains.


%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	depthInsideImage = varargin{iArg}; else
	depthInsideImage = 30;
end, iArg = iArg +1;
if nArg >= iArg && ~isempty(varargin{iArg})
	npts = varargin{iArg}; else
	npts = 350;
end, iArg = iArg + 1;
if nArg >= iArg && ~isempty(varargin{iArg})
    cfg = varargin{iArg};
else
    cfg = struct('debug', false, 'path_data', '.', 'plot_filename', 'domain.png', 'plot_title', '');
end

%
% Reading in "I"
%
if ischar(I)
	I = imread(I);
end
[m,n,o] = size(I);


% perimeter calculation
dx = diff(xp);              dy = diff(yp);  
ds = sqrt(dx.^2 + dy.^2);   perim = sum(ds);   
ss = 2*[0;cumsum(ds)]/perim - 1;
ss(end) = 1;

% find inner boundary
[xp_in,yp_in]  = functions.find_normals(xp,yp,Yhatmax=depthInsideImage);   
%{
[~,~,xc,yc,] = ellipse_fit(xp,yp);
[x_in,y_in] = find_normals(xp,yp,xc,yc);
%}


if cfg.debug
    figure('Visible', 'off');
    imshow(max(I, [], 3), []);
    hold on;
    plot(xp, yp, 'g*');
    plot(xp_in, yp_in, 'r*');
    for i = 1:length(xp)
        plot([xp(i), xp_in(i)], [yp(i), yp_in(i)], 'w-');
    end
    hold off;
    title(cfg.plot_title);
    saveas(gcf, fullfile(cfg.path_data, ['domain_border_', cfg.plot_filename]));
    close(gcf);
end


%
% Building raw data
%
% defining image coordinates
x_im = (1:n);              y_im = (1:m)';
X    = repmat(x_im,m,1);   Y    = repmat(y_im,1,n);


% calculating variables for use in loop
deltax          = diff(xp);
deltay          = diff(yp);
deltax_across   = xp_in(1:end-1) - xp(1:end-1);
deltay_across   = yp_in(1:end-1) - yp(1:end-1);
deltax_across1  = xp_in(2:end) - xp(1:end-1);
deltay_across1  = yp_in(2:end) - yp(1:end-1);


Theta    =  atan2(deltay,deltax);
Xhat1    =  deltax.*cos(Theta) + deltay.*sin(Theta);
Xhat1_in =  deltax_across1.*cos(Theta) + deltay_across1.*sin(Theta);
Yhat1_in = -deltax_across1.*sin(Theta) + deltay_across1.*cos(Theta);
Xhat0_in =  deltax_across.* cos(Theta) + deltay_across.* sin(Theta);
Yhat0_in = -deltax_across.* sin(Theta) + deltay_across.* cos(Theta);
m1       =  Yhat0_in./Xhat0_in;
m2       =  Yhat1_in./(Xhat1_in - Xhat1);
m3       =  (Yhat1_in - Yhat0_in)./(Xhat1_in - Xhat0_in);
c3       =  Yhat0_in - m3.*Xhat0_in;


raw = zeros(npts,o);
for i = 1:npts
    
	% Defining window that we care about 
    % Don't want to do all these calculations on the whole image
	jx = round(xp(i));
	iy = round(yp(i));
	j1 = max(jx-(depthInsideImage+10),1);
	j2 = min(jx+(depthInsideImage+10),n);
	i1 = max(iy-(depthInsideImage+10),1);
	i2 = min(iy+(depthInsideImage+10),m);
	
	% Transforming coordinates on our little local rectangle - image
	theta  = Theta(i);
	Xhat   =  (X(i1:i2,j1:j2) - xp(i))*cos(theta) + (Y(i1:i2,j1:j2) - yp(i))*sin(theta);
	Yhat   = -(X(i1:i2,j1:j2) - xp(i))*sin(theta) + (Y(i1:i2,j1:j2) - yp(i))*cos(theta);
        
	u  = Yhat < m3(i)*Xhat+c3(i) & Yhat >= 0;
	v  = Xhat > (1/m1(i))*Yhat & Xhat < (1/m2(i))*Yhat + Xhat1(i);
	bw = u & v;	
    for j = 1:o
		I1       = I(i1:i2,j1:j2,j);
		I1       = I1(bw);
		raw(i,j) = mean(I1(:));
    end
end


%
% Smoothing. 
%
t = zeros(npts+1,o);
for j = 1:o
	I1      = [raw(1:end-1,j);raw(:,j);raw(2:end,j)];
	smth1   = movmean(I1, round(0.02*npts));           % averaging over ~2% of points
	t(:,j)  = smth1(npts:2*npts);
end
raw = raw([1:end,1],:);
s   = linspace(-1,1,npts+1)';
t   = interp1(ss,t,s);
raw = interp1(ss,raw,s);

if cfg.debug
    figure('Visible', 'off');
    hold on;
    for j = 1:o
        plot(s, t(:,j));
    end
    hold off;
    xlabel('Position (s)');
    ylabel('Intensity');
    title(cfg.plot_title);
    saveas(gcf, fullfile(cfg.path_data, ['domain_signal_', cfg.plot_filename]));
    close(gcf);
end

