function [ DTOut ] = getDTL(Temp, Conc)
%global decaylength %DTOutSave
% global DTliq
    % compute the therman diffusion constant in the liquid
    % IN: Conc(1:K)
    

    K=length(Temp);
    if (length(Conc)~=K )
        error('crap');
    end

KTliq = getKappaTherm(Temp,ones(size(Temp)))*10; % make the diffusion constant consistant with 
                            %the thermal diffusivity it was normalized by
                            %10 so we multiply by 10 to make it correct
                            % dimensions state = 0 give liquid
Cp = 2414e3; %J/(K m^3) % specific heat at melting temp of crystaline
Dbar = 1e-7; %m/s factor by which we normalize diffusion constants

DTliq = KTliq./(Cp.*Dbar);


% DTOut = DTliq*ones(K,1);

DTOut = DTliq;


end

