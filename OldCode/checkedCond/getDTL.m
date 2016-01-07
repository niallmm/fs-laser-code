function [ DTOut ] = getDTL(Temp, Conc)
%global decaylength %DTOutSave
    % compute the therman diffusion constant in the liquid
    % IN: Conc(1:K)
    

    K=length(Temp);
    if (length(Conc)~=K )
        error('crap');
    end

    %DTout = zeros(K,1);
    
 %   DTOut = ones(K,1);

%    DTOut = 5.8e2*ones(K,1);
%    DTOut = Dtherm*DTOut./D0;

    % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
%     
%         DTOut(ii,1) = 1.0;
%     end

% CHECK
% DTOut = 1e3*exp(-Temp/1e0)-1e3.*Conc.^3;
 %============================================
 % Best parameters: D(T) = k_T(T)/(c_p);
 % from a fit of data from Hoglund:
 %  DT = 3.362*exp(-4.5e-3 T) + 0.0824  
 %  DTnorm = DT/Dbar;
 %  Dbar = 1e-7;
 %============================================
  Tmelt = 1685 ; % melting temperature of silicon
  Tambiant = 300  ; % ambient temperature of the substrate
%  
%  T = Temp*(Tmelt-Tambiant)+Tambiant;
%  Dbar = 1e-7;
%  DTOut = (3.362*exp(-T*4.5e-3)*1e-4 + 0.0824e-4)./Dbar;
decaylength = 4.509e-3;
aparam = 3.362e3;
bparam = -decaylength*Tambiant;
cparam = -decaylength*(Tmelt - Tambiant);
dparam = 8.24e1;

DTOut = aparam * exp(bparam + cparam*Temp) + dparam;
%DTOutSave = DTOut;



end

