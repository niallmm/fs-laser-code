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

  
KTsol = getKappaThermdT(Temp,State)*10; % make the diffusion constant consistant with 
                            %the thermal diffusivity it was normalized by
                            %10 so we multiply by 10 to make it correct
                            %dimensions
                            
Cp = 2414e3; %J/(K m^3) % specific heat at melting temp of crystaline
Dbar = 1e-7; %m/s factor by which we normalize diffusion constants

DTdTOut = KTsol./(Cp.*Dbar);
%     DTdTOut = zeros(K,1);



end

