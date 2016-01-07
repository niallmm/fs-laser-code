function [ KappaThermOut ] = getKappaTherm(Temp, State)
    % compute the therman diffusion constant
    % IN: Temp(1:K)
    % IN: Conc(1:K)
    % IN: State(1:K)

    K=length(Temp);
    if (length(State)~=K)
        error
    end

   % KappaThermOut = zeros(K,1);
    KappaThermOut = ones(K,1);
    
    
    
    % ==========================================
    %equilibrium cond
    % ==========================================
%     KappaTherm =1.7e-3 ; % Thermal conductivity
%     KappaThermOut = KappaTherm*KappaThermOut;
%     
    
    
    % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
%     
%     if     (State(ii) ==0) % liquid
%         KappaThermOut(ii,1) = 14.0;
%     elseif (State(ii) ==1) % amorpheous solid
%         KappaThermOut(ii,1) = 0.26;
%     elseif (State(ii) ==2) % crystalline solid
%         KappaThermOut(ii,1) = 2.2;
%     end
% 
%     end

%KappaThermOut = 1.0+1e1*Temp;


% =====================================================================
% Best Parameters: only for crystaline and liquid does not use amorpheous
% k_T = 1/(0.01*(0.033 + 1.55e-3 T + 1.66e-6 T.^2)); in J/(s m K)
% k_T = 4.356*exp(-T*4.18e-3) + 0.2338;  from fit of Hoglund data
% KappaThermOut = k_T/k_Tbar;
% k_Tbar = 10; J/(s m K)
% =====================================================================
Tmelt = 1685 ; % melting temperature of silicon
Tambiant = 300  ; % ambient temperature of the substrate
% 
% T = Temp*(Tmelt-Tambiant)+Tambiant;
% KappaThermOut = 43.56*exp(-T*4.18e-3) + 2.338;

aparam = 43.56;
bparam = -4.181e-3*Tambiant;
cparam = -4.181e-3*(Tmelt-Tambiant);
dparam = 2.338;

KappaThermOut = aparam*exp(bparam+cparam*Temp) + dparam;

end

