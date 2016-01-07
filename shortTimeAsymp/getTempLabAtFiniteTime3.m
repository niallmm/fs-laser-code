function [TempVec,hFinal,Linfty, Umean] = getTempLabAtFiniteTime2( TimeFinal)


global plotflagIC % use that to have the program plot all functions that are integrated numerically
global Dconst_STA alphaParam_STA Linfty_STA maxRelError_STA


%plotflagIC = true;
%plotflagIC = false;


% Set parameters
% Physical parameters
% Dconst = Dconst_STA;%82.4;%1.0; % Heat Diffusion constant 
Dconstl = getDTL(1,0);
Dconsts = getDTS(1,0);
Dconst = [Dconstl Dconsts];
alphaParam = alphaParam_STA;%1.498; % SPECIFY
% Numerical parameters
Linfty = Linfty_STA;%30; % cutoff for the x \to ininity integrals
maxRelError = maxRelError_STA;%1e-8; % Error bound for: abs int_Linfty^{2 Linfty} / int_0^Linfy



% Set time (small) for which we want the approximate solution
% assumption: h(t) = MeanU * t \in [0,TimeFinal]




% Step 1 - solve for the mean velocity Umean

rhs = @(U)   getTermA(U, TimeFinal,Dconst,Linfty,maxRelError) ...
           + getTermB(U,TimeFinal,Dconst,alphaParam) ...
           - getFctJ(U) ;
       
       Uguess = 5;
       Umean = fzero(rhs,Uguess);

       if (Umean <0) 
          Umean 
          error('Umean nagative!!!!')
           
       end
       
if plotflagIC == true
Umean           
tempUplot = linspace(0,3*Umean,10);

% Hacking cuz rhs is not defined as a function operating in vectors - can
% be done but why...
for j=[1:10]
    tempplot(j) = rhs(tempUplot(j));
end
figure(3)       
plot(tempUplot,tempplot);
title('Solving for Umean - find root')
xlabel('Umean');
end

% Step 2 - compute the Temperature at time Tfinal, with h = Umean*TimeFinal

Temperature = @(x) getTermC(x,TimeFinal,Dconst,Linfty,maxRelError) ...
                 + getTermD(x,Umean, TimeFinal, Dconst,alphaParam);

             
% Vectorize that (don't know if that speeds things up but at least it's easier to use)
TempVec = @(xvec) cellfun(Temperature,num2cell(xvec));


% Collect results:

hFinal = Umean*TimeFinal;

% plot temp
  if plotflagIC == true
  

 pltxtemp = linspace(0,3*hFinal,100);
 pltTtemp = TempVec(pltxtemp);
 figure(4)
 plot(pltxtemp,pltTtemp);             
 hold on
 plot([hFinal,hFinal],[0,1.5],'-.');
 hold off
 title('Temperatureprofile');
 xlabel('x');
 ylabel('T');
  end
return

