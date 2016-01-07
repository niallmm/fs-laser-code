function [concLabSurfaceLayer] = getSurfaceLayerGaussianApprox(IntConc,thickness,xlabmesh)


AconcIC = 2 * IntConc ./thickness ./sqrt(pi);

concLabSurfaceLayer = AconcIC*exp(-(xlabmesh.*xlabmesh) ./(thickness).^2);
return
