function [ DTdTOut ] = getDTSdT(Temp, State)
%global decaylength
% ASSUME: The diffusion constant in the SOLID is independent of the solute
% concentration

% compute the derivative of the thermal diffusion constant w.r.t. Temp.
    % IN: Temp(1:K)
    % IN: State(1:K)

  
    K=length(Temp);
%    if (length(Conc)~=K | length(State)~=K)
    if (length(State)~=K)
        error('crap');
    end

 
%    DTdTOut = zeros(K,1);

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

%============================================
 % Best parameters: D(T) = k_T(T)/(c_p);
 % from a fit of data from Hoglund:
 %  DT = (3.362e-4)*exp(-4.5e-3 T) + 8.24e-6 m^2/s 
 %  dDTdT = (1.51e-6)*exp(-4.5e-3 T)
 %  DTnorm = DT/Dbar;
 %  Dbar = 1e-7;
 %============================================
 
 Tmelt = 1685 ; % melting temperature of silicon
 Tambiant = 300  ; % ambient temperature of the substrate
%  
%  T = Temp*(Tmelt-Tambiant)+Tambiant;
%  dTdTemp = Tmelt-Tambiant;
%  Dbar = 1e-7;
%  
% DTdTOut = -(1.51e-6)*exp(-T*4.5e-3)*dTdTemp/Dbar;
decaylength = 4.509e-3;
aparam = 3.362e3;
bparam = -decaylength*Tambiant;
cparam = -decaylength*(Tmelt - Tambiant);
%dparam = 8.24e1;

DTdTOut = aparam * exp(bparam + cparam*Temp) * cparam;

end

