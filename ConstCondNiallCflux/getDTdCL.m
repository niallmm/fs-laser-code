function [ DTdCOut ] = getDTdCL(Temp, Conc)
    % compute the partial derivative of the thermanl diffusion constant in
    % the liquid with respect to concentration
    % IN: Conc(1:K)
    
    %check to make sure Temp and Conc vectors are the same length
    K=length(Temp);
    if (length(Conc)~=K )
        error('crap');
    end

    DTdCOut = zeros(K,1);


end

