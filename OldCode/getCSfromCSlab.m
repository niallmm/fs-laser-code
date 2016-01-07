function [ ConcSolid ] = getCSfromCSlab(h,ConcIntS)
global K1 K2 xmesh concLab xlabmesh K3  concLabbulk Lright
error
% compute the solute concentration in the solid phase at the rescaled grid
% points
% IN: position h

ConcSolid = zeros(K2,1);

% Later: speed up with a fancy vector notation
index = find(h*xmesh(K1+1:K1+K2)> Lright,1,'first');

indexXlargerh  = find(xlabmesh>h,1,'first');

if (isempty(index)) 
%   ConcSolid(1:K2) = ppval(concLabpp,h*xmesh(K1+1:K1+K2));
    ConcSolid(1:K2) = interp1([h;xlabmesh(indexXlargerh:K3)],[ConcIntS;concLab(indexXlargerh:K3)],h*xmesh(K1+1:K1+K2),'spline');
else
    maxindex = index -1;
    %ConcSolid(1:maxindex) = ppval(concLabpp,h*xmesh(K1+1:K1+maxindex));
    ConcSolid(1:maxindex) = interp1([h;xlabmesh(indexXlargerh:K3)],[ConcIntS;concLab(indexXlargerh:K3)],h*xmesh(K1+1:K1+maxindex),'spline');
    ConcSolid(maxindex+1:K2) = concLabbulk;

end


end

