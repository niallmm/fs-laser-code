function [ DTOut ] = getDTS(Temp, State)
%global decaylength
% ASSUME: The diffusion constant in the SOLID is independent of the solute
% concentration
% compute the therman diffusion constant
    % IN: Temp(1:K)
    % IN: State(1:K)
%    global tau0 L0
  
    K=length(Temp);
%    if (length(Conc)~=K | length(State)~=K)
    if (length(State)~=K)
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
%  DTOut = 5.0e2*ones(K,1);

 %============================================
 % Best parameters: D(T) = k_T(T)/(c_p);
 % from a fit of data from Hoglund:
 %  DT = (3.362e-4)*exp(-4.5e-3 T) + 8.24e-6 m^2/s 
 %  DTnorm = DT/Dbar;
 %  Dbar = 1e-7;
 %============================================
 Tmelt = 1685 ; % melting temperature of silicon
 Tambiant = 300  ; % ambient temperature of the substrate
 
%  T = Temp*(Tmelt-Tambiant)+Tambiant;
%  Dbar = 1e-7;
%  
%  DTOut = ((3.362e-4)*exp(-T*4.5e-3) + 8.24e-6)./Dbar;
decaylength = 4.509e-3;
aparam = 3.362e3;
bparam = -decaylength*Tambiant;
cparam = -decaylength*(Tmelt - Tambiant);
dparam = 8.24e1;

DTOut = aparam * exp(bparam + cparam*Temp) + dparam;

end

