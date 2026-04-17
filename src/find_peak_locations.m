function [H,I,bcoor,score] = find_peak_locations(t,nPeaks,s0)

h   = 0.10;

minPeakProminence = 100;
% minPeakDistance   = min(diff(s0));
s = linspace(0,1,length(t));
plot(s,t);
[H,I] = findpeaks(t,'MinPeakProminence',minPeakProminence);
s = linspace(0,1,length(t));
plot(s,t); hold on; plot(s(I),t(I),'*')


while (length(I) < nPeaks)
   minPeakProminence = minPeakProminence - ceil(0.25*minPeakProminence);
   [H,I]             = findpeaks(t,'MinPeakProminence',minPeakProminence); 
end
if (length(I) < nPeaks)
    H = nan; I = nan; bcoor = nan; score = nan;
    disp('Could not required number of peaks in atleast one side of the embryo in the data')
    return
end



%
% Estimate the most probable peak locations
%
% We do this by asking the question - What group of 1:nPeaks peaks are
% most likely the group of nPeaks peaks for this gene? We develop a
% score based on relative peak location, and how far peaks are from
% canonical peak locations 
%
extraPeaks  = length(I) - nPeaks;
[I,i0]      = sort(I);         
H           = H(i0);

ds0         = zeros(extraPeaks+1,1);
dlocs       = ds0;
for i=1:extraPeaks+1
    locs     = s(I(i:i+nPeaks-1));   
    dlocs(i) = var(locs);
    ds0(i)   = norm(s0 - locs);    
end 
[~,i0]  = min(ds0+dlocs);
H       = H(i0:nPeaks+i0-1);
I       = I(i0:nPeaks+i0-1);
score   = norm(s0-s(I));



% find initial estimates of the boundaries
%
bcoor   = zeros(nPeaks,2); 
idx     = find((t-h*H(1)).*(t([2:end,1])-h*H(1)) < 0);
jf      = find(idx < I(1));
if isempty(jf)
    jf = 1;
end
jf             = jf(end);
loc            = 0.5*(s(idx(jf))   + s(idx(jf)+1));
[~,bcoor(1,1)] = min(abs(s-loc));


% check by plotting
% This is the first point!
%{
figure
plot(s,(t-h*H(1)).*(t([2:end,1])-h*H(1)))
hold on
plot(s(bcoor(1,1)),t(bcoor(1,1)),'*')
%}


if nPeaks>1
    for i = 2:nPeaks
        t2              = t(I(i-1):I(i));
        [~,i0]          = min(t2);
        bcoor(i-1,2)    = I(i-1) + i0;
        bcoor(i,1)      = bcoor(i-1,2);         
     end
end

idx = find((t-h*H(end)).*(t([2:end,1])-h*H(end)) < 0);
jf  = find(idx < I(end));
if isempty(jf)
    jf = 1;
end
jf = jf(end);
if jf == length(idx) || isempty(idx)
    bcoor(end,2) = length(s);
else
    if length(s) > idx(jf+1)+1
        loc         = 0.5*(s(idx(jf+1)) + s(idx(jf+1)+1));
        [~,bcoor(end,2)]= min(abs(s-loc));
    else
        bcoor(end,2) = length(s);
    end
end











































end