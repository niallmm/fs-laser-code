function [ KappaThermdCOut ] = getKappaThermdC(Temp, State)
    % compute the der. of therman diffusion constant wrt conc
    % IN: Temp(1:K)
    % IN: Conc(1:K)
    % IN: State(1:K)

    K=length(Temp);
    if (length(State)~=K)
        error
    end

    KappaThermdCOut = zeros(K,1);


end

