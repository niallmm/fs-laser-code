function [ KappaThermdTOut ] = getKappaThermdT(Temp, State)
    % compute the der. of therman diffusion constant wrt temp
    % IN: Temp(1:K)
    % IN: Conc(1:K)
    % IN: State(1:K)

    K=length(Temp);
    if (length(State)~=K)
        error
    end

%    KappaThermdTOut = zeros(K,1);

    % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
%     
%     if     (State(ii) ==0) % liquid
%         KappaThermdTOut(ii,1) = 0.0;
%     elseif (State(ii) ==1) % amorpheous solid
%         KappaThermdTOut(ii,1) = 0.0;
%     elseif (State(ii) ==2) % crystalline solid
%         KappaThermdTOut(ii,1) = 0.0;
%     end
% 
%     end


% KappaThermdTOut = 1e1;


% =====================================================================
% Best Parameters: only for crystaline and liquid does not use amorpheous
% k_T = 1/(0.01*(0.033 + 1.55e-3 T + 1.66e-6 T.^2)); in J/(s m K)
% k_T = 4.356*exp(-T*4.18e-3) + 0.2338;  from fit of Hoglund data
% KappaThermOut = k_T/k_Tbar;
% k_Tbar = 10; J/(s m K)
% =====================================================================
Tmelt = 1685 ; % melting temperature of silicon
Tambiant = 300  ; % ambient temperature of the substrate
% dTdTemp = Tmelt-Tambiant;
% T = Temp*(Tmelt-Tambiant)+Tambiant;
% KappaThermdTOut = -0.182*exp(-T*4.18e-3)*dTdTemp;


aparam = 43.56;
bparam = -4.181e-3*Tambiant;
cparam = -4.181e-3*(Tmelt-Tambiant);
%dparam = 2.338;

KappaThermdTOut = aparam*exp(bparam+cparam*Temp) * cparam;

% end

