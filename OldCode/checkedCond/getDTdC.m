function [ DTdCOut ] = getDTdC(Temp, Conc, State)
error
% compute the derivative of the thermal diffusion constant w.r.t. Temp.
    % IN: Temp(1:K)
    % IN: Conc(1:K)
    % IN: State(1:K)

  
    K=length(Temp);
    if (length(Conc)~=K | length(State)~=K)
        error('crap');
    end

    DTdCOut = zeros(K,1);

    % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
%     
%     if     (State(ii) ==0) % liquid
%         DTdCOut(ii,1) = 0.0;
%     elseif (State(ii) ==1) % amorpheous solid
%         DTdCOut(ii,1) = 0.0;
%     elseif (State(ii) ==2) % crystalline solid
%         DTdCOut(ii,1) = 0.0;
%     end
% 
%     end
%DTdCOut = 2.e1*Conc;% + zeros(K,1);

end
