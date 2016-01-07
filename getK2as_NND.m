function [K2as] = getK2as_NND(oldh,tglobal)
global K1 K2 Dconst_BBC xmesh TintInvLam_global
largeparam = 10;
h_threshold = 100;
K2as = find(oldh*xmesh(K1+1:K1+K2)- largeparam * sqrt(4*Dconst_BBC*tglobal)>0,1,'first');

if(~isempty(K2as))
    one_over_lambda = TintInvLam_global(oldh*xmesh(K2as));
    t_threshold =1./largeparam ./ one_over_lambda ./one_over_lambda ./Dconst_BBC;
    
    if((tglobal>t_threshold)|| (oldh > h_threshold))
       K2as = K2; 
    end
        
    
    
% end
% if (isempty(K2as) || (oldh>2))
else
    K2as = K2 ;
end
%K2as = K2-1;
%K2as
return
