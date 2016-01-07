function [ DCdCOut ] = getDCdCL(Temp, Conc)
    % compute the der. of solute diffusion constant wrt to conc (liquid)
    % IN: Temp(1:K)
    % IN: Conc(1:K)

    K=length(Temp);
    if (length(Conc)~=K)
        error
    end

    DCdCOut = zeros(K,1);

    % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
% 
%         DCdCout(ii,1) = 0.0;
% 
%     end
%DCdCOut = 15e1*Conc.^2;

end

