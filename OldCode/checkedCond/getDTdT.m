function [ DTdTOut ] = getDTdT(Temp, Conc, State)
error
% compute the derivative of the thermal diffusion constant w.r.t. Temp.
    % IN: Temp(1:K)
    % IN: Conc(1:K)
    % IN: State(1:K)

  
    K=length(Temp);
    if (length(Conc)~=K | length(State)~=K)
        error('crap');
    end

 
    DTdTOut = zeros(K,1);

%     % Speed that up with a clever vector notation ...
%     
%     for ii = 1:K
%     
%     if     (State(ii) ==0) % liquid
%         DTdTOut(ii,1) = 0.0;
%     elseif (State(ii) ==1) % amorpheous solid
%         DTdTOut(ii,1) = 0.0;
%     elseif (State(ii) ==2) % crystalline solid
%         DTdTOut(ii,1) = 0.0;
%     else
%         error('state unknown in getDTdT')
%     end
% 
%     end
%DTdTOut = 1e3+zeros(K,1);

end

