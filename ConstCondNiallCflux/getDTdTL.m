function [ DTdTOut ] = getDTdTL(Temp, Conc)
% global decaylength
% compute the partial derivative of the thermanl diffusion constant in
    % the liquid with respect to Temperature
    % IN: Conc(1:K)
    

    K=length(Temp);
    if (length(Conc)~=K )
        error('crap');
    end

    DTdTOut = zeros(K,1);


end

