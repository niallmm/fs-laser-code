function [ KappaThermdTOut ] = getKappaThermdT(Temp, State)
    % compute the der. of therman diffusion constant wrt temp
    % IN: Temp(1:K)
    % IN: Conc(1:K)
    % IN: State(1:K)

    K=length(Temp);
    if (length(State)~=K)
        error
    end
    
    Tempvec = [273 300 350 400 500 600 700 800 900 1000 1100 1200 1300 1400 ...
        1500 1600 1685];
    Kappavec = [1.68 1.48 1.19 0.989 0.762 0.619 0.508 0.422 0.359 0.312 ...
        0.279 0.257 0.244 0.235 0.227 0.221 0.220]*10;
    delta = [diff(Kappavec)./diff(Tempvec) 0];
    Temp2 = [(Tempvec(1:end-1)+Tempvec(2:end))/2 2500];

    KappaThermdTOut = interp1(Temp2, delta, Temp, 'linear', 'extrap');


end

