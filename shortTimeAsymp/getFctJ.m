function [temp] = getFctJ(vel)
global ParamMu_STA
ParamMu = ParamMu_STA;%37.0;

temp = 1.0 + vel./ParamMu;

return