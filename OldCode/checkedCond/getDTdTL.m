function [ DTdTOut ] = getDTdTL(Temp, Conc)
% global decaylength
% compute the partial derivative of the thermanl diffusion constant in
    % the liquid with respect to Temperature
    % IN: Conc(1:K)
    

    K=length(Temp);
    if (length(Conc)~=K )
        error('crap');
    end

%    DTdTOut = zeros(K,1);

    % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
%     
%         DTdTOut(ii,1) = 0.0;
%     end

%     DTdTOut = -1e3*exp(-Temp/1e0);
%============================================
 % Best parameters: D(T) = k_T(T)/(c_p);
 % from a fit of data from Hoglund:
 %  DT = 3.362*exp(-4.5e-3 T) + 0.0824  
 %  dDTdT = 0.0151*exp(-4.5e-3 T)
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
%  
% DTdTOut = 0.0151*exp(-T*4.5e-3)*dTdTemp*1e-4/Dbar;
% 
decaylength = 4.509e-3;
aparam = 3.362e3;
bparam = -decaylength*Tambiant;
cparam = -decaylength*(Tmelt - Tambiant);
%dparam = 8.24e1;

DTdTOut = aparam * exp(bparam + cparam*Temp) * cparam;

end

