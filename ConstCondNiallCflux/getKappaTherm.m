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
% 
Tmelt = 1685;    %K
% Tmelt = 1430; %K amorphous
Tambiant = 300;  %K   
% Tempvec = [273 300 350 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1685];
% Kappavec = [1.68 1.48 1.19 0.989 0.762 0.619 0.508 0.422 0.359 0.312 0.279 0.257 0.244 0.235 0.227 0.221 0.220]*10;
 a = 43.56;
b = -0.004181;
c = 2.338;
% Kout  = a*exp(b*Tempn)+c;
KappaThermOut(State==1) = 14;
dimenT = (Temp(State==2)*(Tmelt-Tambiant)-Tambiant);
KappaThermOut(State==2) = a*exp(b*dimenT)+c;


        


end

