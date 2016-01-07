function [ DTOut ] = getDT(Temp, Conc, State)
error
% compute the therman diffusion constant
    % IN: Temp(1:K)
    % IN: Conc(1:K)
    % IN: State(1:K)
%    global tau0 L0
  
    K=length(Temp);
    if (length(Conc)~=K | length(State)~=K)
        error('crap');
    end
 
    %DTOut = zeros(K,1);
    
    
    % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
%     
%     if     (State(ii) ==0) % liquid
%         DTOut(ii,1) = 1.0;
%     elseif (State(ii) ==1) % amorpheous solid
%         DTOut(ii,1) = 1.1;
%     elseif (State(ii) ==2) % crystalline solid
%         DTOut(ii,1) = 1.2;
%     end
% 
%     end
% CHECK
%DTOut = 1+1e1*Conc.^2 + 1e3*Temp;%+ exp(Conc);

% NOTE - HAS TO DEPEND ON TERPERATURE
 DTOut = 5.0e2*ones(K,1);
  

end

