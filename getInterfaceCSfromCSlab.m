function [ ConcInterfaceSolid ] = getInterfaceCSfromCSlab(h)
global K1 K2 L1 xmesh concLabpp concLabbulk Lright
% compute the solute concentration in the solid phase at the rescaled grid
% points
% IN: 



if (Lright>h*L1) 
   ConcInterfaceSolid = ppval(concLabpp,h*L1);
else
   ConcInterfaceSolid = concLabbulk;

end


end

