function [xlabmesh] = getXlabMesh(K3,Lright,K1,hin)

dx = (Lright-hin)./(K3-K1);

xlabmesh = [hin*linspace(0,1.0,K1)';linspace(hin+dx,Lright,K3-K1)'];

return