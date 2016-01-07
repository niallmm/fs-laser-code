function [ DTdCOut ] = getDTdCL(Temp, Conc)
    % compute the partial derivative of the thermanl diffusion constant in
    % the liquid with respect to concentration
    % IN: Conc(1:K)
    

    K=length(Temp);
    if (length(Conc)~=K )
        error('crap');
    end

    DTdCOut = zeros(K,1);

    % Speed that up with a clever vector notation ...
    
%     for ii = 1:K
%     
%         DTdCOut(ii,1) = 0.0;
%     end

%    DTdCOut = zeros(K,1) -3e3.*Conc.^2;

end

