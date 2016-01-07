function [xmesh] = GetSolidGridVers2(xinf,Npoints,param1, param2, dximpose)

temp = linspace(0,1,Npoints);

L2 = xinf;
bparam = param1;%5.0; % impose that
cparam = param2;%4.0; % impose that

slope = dximpose * Npoints; % impose that





 rhs = @(d) L2 - slope./d * exp(bparam+d) -1.0 + slope./(d);
 
 dguess = 1.0;
 
 dzero = fzero(rhs,dguess);
 
 dparam = dzero;
 aparam = slope ./ dparam;
 eparam = 1.0 - aparam;
 
 
 xmesh = aparam .* exp(bparam .* temp.^cparam + dparam .*temp) + eparam;
 xmesh = xmesh';
 size(xmesh)
 


return