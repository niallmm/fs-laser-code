function [ StateInterfaceSolid ] = getInterfaceStateFromstateLab(h)
global K1 L1 K2 xmesh stateLab xlabmesh stateLabbulk
% compute the solute concentration in the solid phase at the rescaled grid
% points
% IN: 


%StateInterfaceSolid = interp1([xlabmesh;1.001.*xlabmesh(end);1e10],[stateLab;stateLabbulk;stateLabbulk],h*L1,'nearest');
StateInterfaceSolid = myinterp1nearest([xlabmesh;1.001.*xlabmesh(end);1e10],[stateLab;stateLabbulk;stateLabbulk],h*L1);

end
