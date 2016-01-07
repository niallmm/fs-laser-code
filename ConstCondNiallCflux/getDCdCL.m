function [ DCdCOut ] = getDCdCL(Temp, Conc)
    % compute the der. of solute diffusion constant wrt to conc (liquid)
    % IN: Temp(1:K)
    % IN: Conc(1:K)

    K=length(Temp);
    if (length(Conc)~=K)
        error
    end

    DCdCOut = zeros(K,1);


end

