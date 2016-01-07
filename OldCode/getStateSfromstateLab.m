function [ StateSolid ] = getStateSfromstateLab(h)
global K1 K2 xmesh stateLab xlabmesh stateLabbulk
% compute the solute concentration in the solid phase at the rescaled grid
% points
% IN: 
if (h<0) 
    error('h< 0 in getStateSfromstateLab')
end
StateSolid = zeros(K2,1);

%StateSolid(1:K2) = interp1([xlabmesh;1.001.*xlabmesh(end);1e10],[stateLab;stateLabbulk;stateLabbulk],h*xmesh(K1+1:K1+K2),'nearest');
StateSolid(1:K2) = myinterp1nearest([xlabmesh;1.001.*xlabmesh(end);1e10],[stateLab;stateLabbulk;stateLabbulk],h*xmesh(K1+1:K1+K2));

end
