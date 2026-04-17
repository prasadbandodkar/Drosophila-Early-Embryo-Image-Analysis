function [cfun,cvals,cint,cint68,r2] = fitelephant(t,s,A,I,f)

k    = isnan(A);
A(k) = min(A(~k));
I(k) = 1;

%
% Upper and lower limits for the parameters (except x0)
%
A(A == 0) = 1e-2;
AL = 0.01*A;        AU = 10*A;
BL = 0;             BU = min(A);        B = 0.5*BU;     % inital bckgrnd is half of min(A)
delt  = ones(1,length(A)); 
deltL = 0.5*delt; 
deltU = 2*delt;

%
% Determining what the values for x0 (initial guesses as well as upper and
% lower bounds) should be.  These change depending on whether we have
% lateral genes or not.  Lateral genes will have some intermediate value
% for x0 and will have upper and lower bounds.  Fully dorsal genes and
% fully ventral genes do not have an x0 for fitting purposes.
%
cnames = coeffnames(f);
nx     = sum(functions.strfindDU(cnames,'x'));
if nx == length(I)
	x0 = s(I);
else
	x0 = zeros(nx,1);
	count = 1;
	for i = 1:length(I)
		if I(i) ~=1 && I(i) ~= length(s) && ~k(i)
			x0(count) = s(I(i));
			count = count + 1;
		end
	end
end
x0L = x0 - 0.2; x0U = x0 + 0.2; % Searching +-20% from x0 (should this come from sB?)

%
% Implementing our options
%

opts = fitoptions('Method','NonlinearLeastSquares',...
    'Startpoint',[A, B, delt, x0'],...
    'Lower',     [AL, BL, deltL, x0L'],...
    'Upper',     [AU, BU, deltU, x0U']);

%
% The actual fit, followed by unpacking the coefficient values, confidence
% intervals, and r2 value.
%
[cfun,gof] = fit(s,t,f,opts);
cvals      = coeffvalues(cfun);
cint       = confint(cfun);
cint68     = confint(cfun,(1-2*(1-normcdf(1,0,1))));
r2         = gof.rsquare;

%{
filenameshort = cell2mat(c(end-2:end-1));

idx = num2strDU(i,3);
if ~exist('Fittedpeaksimages','dir')
	mkdir Fittedpeaksimages
end
print(456,'-djpeg',['Fittedpeaksimages',filesep,filenameshort,'_',...
	genenames,'_',genotype,idx,'.jpg'])
close(456)
%}

end