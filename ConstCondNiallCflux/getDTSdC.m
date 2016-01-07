function [ DTdCOut ] = getDTSdC(Temp, State)
error('not needed!!!')
% ASSUME: The diffusion constant in the SOLID is independent of the solute
% concentration

% compute the derivative of the thermal diffusion constant w.r.t. Temp.
    % IN: Temp(1:K)
    % IN: State(1:K)

  
    K=length(Temp);
    if (length(State)~=K)
%    if (length(Conc)~=K | length(State)~=K)
        error('crap');
    end

    DTdCOut = zeros(K,1);



end

