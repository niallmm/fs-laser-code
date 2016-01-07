function [ KappaThermOut ] = getKappaTherm(Temp, State)
    % compute the therman diffusion constant
    % IN: Temp(1:K)
    % IN: Conc(1:K)
    % IN: State(1:K)

    K=length(Temp);
    if (length(State)~=K)
        error
    end

    KappaThermOut = zeros(K,1);

    


% =====================================================================
% From Hoglund:
% k_T = 1.4 J/(K s cm) = 140 J/(K s m) in Liquid
% k_T = 0.22 J/(K s cm) = 22 J/(K s m) in crystaline at crystaline melt temp
% we actually have data for crystaline as a function of temperature, but
% for most of what is happening it is well approximated at melting temp
% KappaThermOut = k_T/k_Tbar;
% k_Tbar = 10; J/(s m K)
% =====================================================================


for ii = 1:K
    if State(ii) ==0
        KappaThermOut(ii) = 14.0; % liquid
    else
        KappaThermOut(ii) =  2.2; % solid crystaline
    end
end
        


end

