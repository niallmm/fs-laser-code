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

    % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
%     
%     if     (State(ii) ==0) % liquid
%         KappaThermdCOut(ii,1) = 0.0;
%     elseif (State(ii) ==1) % amorpheous solid
%         KappaThermdCOut(ii,1) = 0.0;
%     elseif (State(ii) ==2) % crystalline solid
%         KappaThermdCOut(ii,1) = 0.0;
%     end
% 
%     end

end

