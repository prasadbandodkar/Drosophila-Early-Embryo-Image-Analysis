function t = get_side(t,iPoles,side)

arguments
    t      (:,1) double
    iPoles (1,2) double
    side   (1,1) double {mustBeMember(side, [1,2])}
end

ns    = length(t)-1;
apole = round(iPoles(1));
ppole = round(iPoles(2));

t = t(1:ns);
if side == 1        % anterior to posterior or dorsal to ventral
    if apole < ppole
        t = t(apole+1:ppole);
    else
        t = [t(apole+1:ns);t(1:ppole)];
    end
elseif side == 2    % posterior to anterior or ventral to dorsal flipped    
    if apole < ppole
        t = [flipud(t(1:apole-1)); flipud(t(ppole:ns))];
    else
        t = flipud(t(ppole:apole-1));
    end
end


end