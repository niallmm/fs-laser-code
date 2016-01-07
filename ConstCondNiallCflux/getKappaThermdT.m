function [ KappaThermdTOut ] = getKappaThermdT(Temp, State)
    % compute the der. of therman diffusion constant wrt temp
    % IN: Temp(1:K)
    % IN: Conc(1:K)
    % IN: State(1:K)

    K=length(Temp);
    if (length(State)~=K)
        error
    end

    KappaThermdTOut = zeros(K,1);

end

