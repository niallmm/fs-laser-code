function [xmesh] = GetSolidGrid(xinf,Npoints,param)

% temp = linspace(0,1,Npoints);
% B = param;
% expB = exp(B);
% C = (expB - xinf)/(expB-1);
% A=1-C;
% 
% xmesh = A*exp(B.*temp) +C;

% Play around with another version to alow for more resolution close to the
% enforced BC at xinf

Nend = 10;
temp = linspace(0,1,Npoints-Nend);
B = param;
expB = exp(B);
C = (expB - 0.8*xinf)/(expB-1);
A=1-C;

xmeshexpPart = A*exp(B.*temp) +C; % last point is at 0.95 xinf
xmeshend = linspace(0.8*xinf,xinf,Nend+1);
xmesh = [xmeshexpPart,xmeshend(2:end)];

return



% temp = linspace(0,1,Npoints);
% B = param;
% expB = exp(B);
% C = (expB - xinf)/(expB-1);
% A=1-C;
% 
% xmesh = A*exp(B.*temp) +C;
