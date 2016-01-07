function [ DCdTout ] = getDCdTL(Temp, Conc)
    % compute the der. of solute diffusion constant wrt to temp (liquid)
    % IN: Temp(1:K)
    % IN: Conc(1:K)

    % check that the Temp and Conc vectors are the same length
    K=length(Temp);
    if (length(Conc)~=K)
        error
    end



DCdTout = zeros(K,1);


end

