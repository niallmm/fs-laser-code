function [concLabNew,stateLabNew,tend,hverst,StateTime,Time,Tsaved, xliquid,Concliqsave] = OneShotTestAdv(hIC,TIC, concLabIC, concLabbulkIC,stateLabIC, stateLabbulkIC,tIC,OutputTime)

global K1 K2 Z xmesh plotflag
global dt dtmax dhmax convNewton convthreshld maxIterNewton maxSteps deltaTsave min_h_resold
global xlabmesh 
global solidStateGlobal% concSolidEnforceGlobal
global concLab concLabpp stateLab tglobal concLabbulk stateLabbulk 

global K2as_NND time_NND 
fakecounter = 0;

outputTestflag = true;


concLab = concLabIC;
concLabbulk = concLabbulkIC;
stateLab = stateLabIC;
stateLabbulk = stateLabbulkIC;

solidStateIC = getInterfaceStateFromstateLab(hIC); % at the interface

% Build Initial state vector
% 

% generate spline representation of the concentration
concLabpp = spline([xlabmesh],[concLabIC]); % one can prescribe the end slopes!!
% evaluate the spline at the hIC*xmesh(1:K1+1) (including the point in the
% solid
CIC = zeros(K1+1,1);
CIC(1:K1+1) = ppval(concLabpp,hIC*xmesh(1:K1+1));
% Initial state vector
% T[1:K1] liquid, T[K1+1,K1+K2] solid, C[1:K1] liquid, C[K1+1] (SINGLE
% POINT) solid, h
uIC = [TIC;CIC;hIC];






t=tIC;

% Generate the newton iterator

IFunc = crank_newton_generator(@getresid, @getjacob,convNewton,maxIterNewton) ;



% Now the main loop in the forward time integration


h=hIC;
%hbreak = 0.5*hIC;  % value of h where we end the computation
%hbreak = 0.99*hIC;  % value of h where we end the computation
hbreak = min_h_resold;
hbreakforce = 0.99*hIC;
% This is a globally defined variable that contains the (interpolated)
% concentration in the solid at position x=h (enforced during the melting
% process but NOT during resolidification)
% concSolidEnforceGlobal = CIC(K1+1);
% Global variable for the solid state AT the interface
solidStateGlobal = solidStateIC;
% add time as a globally available quantity (for explicitly time dependent
% conditions)
tglobal =t;
K2as_NND = getK2as_NND(hIC,tglobal);


uold=uIC;

solidflag = false; % indicate that we are resolidifying (hdot <0)

size(uIC)
% for output purposes
%hverst = zeros(500,2); % preaccolacte memory
%Tathsaved = zeros(500,3);
tsaveout = t;
savecounter = 1;
plotcounter = 1;
plott = t;

for ii = 1:maxSteps % MAIN LOOP
     %   ii, t,
     % variable timestep version similar to Scotts code
        % first we need to take an accurate step
    accurate = false ;      
    
    
    
    while (~accurate)
        % take two dt/2-steps, then a dt-step
        success = true ;
        if (success) ; 
            time_NND = t+dt/2;
                       size(uold)
            [guess2, success] = IFunc(uold, dt/2) ; 
 
        end ;
        if (success) ; 
            time_NND = t+dt;
            [guess2, success] = IFunc(guess2, dt/2) ;
        end ;
        if (success) ;
            time_NND = t+dt;
            [guess1, success] = IFunc(uold, dt, guess2) ; 
        end ;

        % if successful, then compare for error
        if (success) ;
           
            % CHECK NORM USED FOR ERROR ESTIMATES - IT"S A LOCAL ONE!!!!
          %  err = max( (guess2 - guess1)./(guess2+convthreshld));
          %  err = max( (guess2 - guess1));
          %  err = max( (guess2 - guess1));
            err = max(abs( (guess2(1:K1+K2) - guess1(1:K1+K2))));
%            err = max(abs( (guess2(1:K1+K2as_NND-2) - guess1(1:K1+K2as_NND-2))));
            accurate = (err < convthreshld); %*dt/tfinal) 
  if (plotflag ==1)
    % For testing purposes - plot stuff
    figure(1)       
    subplot(4,1,4)
    plot(abs(guess2-guess1),'k'); % concentration in the solid
    title('error')
    drawnow;    
  
    
    
  end
            % Initical Temperature profile has jumps... maybe reduce the
            % accuracy requirements
        %    if (ii <5) ; accurate = (err < 2e-1) ; end;
        %    if (ii <15) ; accurate = (err < 2) ; end;
        end ;

        if (~success)
                
                fprintf('Newton not converging\n')
        end
        % cut timestep if unsuccessful or inaccurate
        if (~accurate)
            dt = dt / 2 ;
            if (dt < 1e-16) ;  
                fprintf('\n\n   *** %s ***   \n\n', 'dt too small'); return ; 
                if (~success)
                
                fprintf('Newton not converging\n')
                end
            end ;
            
        end 

          

    end ;
    
    % increase dt if we're too accurate
    if (err < 0.02*convthreshld)
        deltah = abs(guess2(end)-h);
        if (deltah <dhmax)
        if (dt < dtmax)
     %   dt = dt * 2^(1/4) ;
        dt = dt * 2^(1/8) ;
        end
        end
    end
    
     % check if the interface moves back (resolidification)
    if (~solidflag)
       if (guess2(end) < h )
       solidflag = true;
       end
    else   
    if (guess2(end) > h )
       solidflag = false;
       end
       
    end
   
    
    h = guess2(end);
    solidconc = guess2(end-1);
   
   
    
    
    
    
    
    
    
    
    
    % update global variables
    % The global variables contain
    % 1. the state (amorpheous / crystalline ) of the solid at the moving
    % front
    % 2. the time: tglobal (for time dependent BC or something like that)
    % 3. the solute concentration at the interface that is enforeced for
    % hdot > 0
    
  
    
    
    if (~solidflag) % melting
        % Concept: hdot > 0 get info from the lab-frame data
        % for hdot > 0 both the conc. and the state are simply given by the lab
        % frame date
        solidStateGlobal = getInterfaceStateFromstateLab(h);
%         concSolidEnforceGlobal = getInterfaceCSfromCSlab(h);
    
    else
        % for hdot < 0 the conc. in the solid is dynamically chosen (see
        % guess2(Z+K1+1)   )
        % and the state in which the solid forms has to be determined by some
        % sort of criterion which we don't know yet.
    
        solidStateGlobal = 2; % LATER: put a fancy external function tha depends on
        % hdot, interface temp, the two concentrations
        
        % update lab frame data
        indexXsmallerh = find(xlabmesh<h,1,'last');
        indexXlargerh  = find(xlabmesh>h,1,'first');
                oldh = uold(end);
        indexXsmallerhold = find(xlabmesh<oldh,1,'last');
        
         if (h < .1)
            solidconc
            
        concLab(indexXlargerh:indexXsmallerhold)  = ...
        interp1([h,xlabmesh(indexXsmallerhold+1)'], ...
                 [solidconc,concLab(indexXsmallerhold+1)'], ...
                 xlabmesh(indexXlargerh:indexXsmallerhold),'linear','extrap');
        
        
        
        % set all grid points in the lab frame < h to the concentration at x=h 
        % or interpolate after moving past the point !
        %concLab(1:indexXsmallerh) = solidconc;
        % for the point we have just passed we interpolate:
        concLab(1:indexXsmallerh) = ...
        interp1([h,xlabmesh(indexXlargerh:indexXlargerh+1)'], ...
                 [solidconc,concLab(indexXlargerh:indexXlargerh+1)'], ...
                 xlabmesh(indexXsmallerh),'linear','extrap');
         else
            
        % set all grid points in the lab frame < h to the concentration at x=h 
        % or interpolate after moving past the point !
        concLab(1:indexXsmallerh) = solidconc;
        
         end
        % for the point we have just passed we interpolate:
     %   concLab(indexXlargerh) = ...
      %      interp1([h,xlabmesh(indexXlargerh+1)],[solidconc,concLab(indexXlargerh+1)],xlabmesh(indexXlargerh));
        
        
        
        % check that
        if( (xlabmesh(indexXlargerh) - xlabmesh(indexXsmallerh)) / ...
                h                    - xlabmesh(indexXsmallerh)) < 0.5
            stateLab(indexXlargerh) = solidStateGlobal;
            
        end
        
    end
    
    
    
    % PLOT AND STORE STUFF
    if (plotflag >0)
    % For testing purposes - plot stuff
    figure(1)       
    subplot(4,1,1)
 %       plot(xmesh(1:K1+K2)*h,guess2(1:K1+K2),'r'); % Temperature
 %       plot(xmesh(1:K1)*h,guess2(1:K1),'r'); % Temperature
 %   plot(xmesh(1:K1+ceil(2*K2/4))*h,guess2(1:K1+ceil(2*K2/4)),'r'); % Temperature
    plot(xmesh(1:K1+K2)*h,guess2(1:K1+K2),'r'); % Temperature
 
     axis([0 10 0 2.5]) 
    title('Temp');
    %hold on;
    subplot(4,1,2)
    plot(xmesh(1:K1+1)*h,guess2(Z+1:Z+K1+1),'b'); % concentrtion in liquid

    title('concentration (liquid)');
    subplot(4,1,3)
    plot(xlabmesh,concLab,'b'); % concentration in the solid
    title('concentration (solid)')
    drawnow;    
  
%     figure(99)
%     plot(guess2(1:K1+K2),'r'); % Temperature
%   
%    min(guess2(1:K1+ceil(2*K2/4)))
    end
    if (plotflag == -1)
           if (t+dt>plott+ 0.0005)
                plott = t+dt;
                solidindex = find(xlabmesh>h);
        videofile = ['video/' sprintf('%03d', plotcounter)];
        figure(1)
        subplot(2,1,1)
        plot(xmesh(1:K1+K2)*h*10,guess2(1:K1+K2)*1385+300,'b','LineWidth',2); % Temperature
        hold on
%         title('Temperature Profile')
        line([h*10; h*10],[1000; 3500], 'LineWidth',2,'Color', 'k','LineStyle', '--')       
        axis([0 40 1000 3500])
%         xlabel('depth [nm]')
        ylabel('Temperature [K]')
        hold off 
        subplot(2,1,2)
        plot(xmesh(1:K1+1)*h*10, guess2(Z+1:Z+K1+1), 'k', 'LineWidth',2)
        hold on
%         title('Concentration Profile')
        xlabel('depth [nm]')
        ylabel('Normalized Concentration')
        if t > 0.1
        plot(xlabmesh(solidindex)*10, concLab(solidindex), 'k', 'LineWidth',2)
        end
        hold off
        line([h*10; h*10],[0; 3], 'LineWidth',2,'Color', 'k','LineStyle', '--')
        axis([0 40 0 3])
        set(gcf,'color','w')
        saveas(gcf, videofile, 'jpg')
        
        plotcounter = plotcounter+1;
           end
    end
        
   if (t+dt>tsaveout+ deltaTsave)
        hverst(savecounter,1) = t+dt;
        hverst(savecounter,2) = guess2(end);
%        Tathsaved(savecounter, 1) = t+dt;
%         Tathsaved(savecounter,2) = guess2(K1);
%         Tathsaved(savecounter,3) = guess2(K1+1);
        Tsaved(savecounter, 1) = t+dt;
        Tsaved(savecounter, 2:K1+K2+1) = guess2(1:K1+K2);
        tsaveout = t+dt;
        xliquid(savecounter,1:K1+K2) = xmesh(1:K1+K2)*h;
        Concliqsave(savecounter,1:K1+K2) = guess2(1:K1+K2);
        
%         global L2
%         fname = sprintf('Temperature_New_K2%d_L2%d_no%d.mat', K2, L2, savecounter)
%         save(fname,'xmesh','h','guess2','tsaveout','K1','K2','Z');
         savecounter=savecounter +1;
     end
    
    
    
    
    
    
    if guess2(end) <hbreakforce
        hbreakforce1 = hbreakforce
        break
    end
    if (solidflag)
        if guess2(end) <hbreak  
            hbreak1 = hbreak
            break
        end
         
    end
   
    
    %
% %     % CHANGE THAT
%      if (solidflag)
%          fakecounter = fakecounter + 1;
% %         
%      end
% %     
%      if (fakecounter >10)
%         break 
%      end
%     
    uold = guess2;
    t = t+dt;
    tglobal = t;
    if(~solidflag)
        K2as_NND = getK2as_NND(guess2(end),tglobal);
    end
    

   if (t>OutputTime && outputTestflag == true)
       Time=t;
       StateTime = uold;
       outputTestflag = false;
      % break
   end
    
    
end
ii

% Return values

concLabNew = concLab;
% Special treatment for the points smaller than hbreak:
     indexXsmallerh = find(xlabmesh<h,1,'last');
     indexXlargerh  = find(xlabmesh>h,1,'first');
     % NEW
    % Interpolate concentration in the solid to h -> 0 based on known
    % asyptotic behavior
    
    h0as = h;
    C0as = solidconc;
    C0liquid = guess2(end-2);
    % extract k(v) from the last two concentrations
    kvas = C0as./C0liquid;
    
    
    % option 1: Real asymptotics BUT divergence

    if(false)
    concLabNew(2:indexXsmallerh) = ...
        C0as .* (h0as./xlabmesh(2:indexXsmallerh)) .^(1-kvas);
    
    concLabNew(1) = concLabNew(2); % cut divergence
    
     
    else  
    % Option 2 - Gaussian approximation satisfying solute conservation
    Has = 0.5 * h0as;
    
    Aas = (1.0 - kvas).* h0as * C0as ./ ...
            (   kvas .* Has .* sqrt(pi) ./2.0 .* erf(h0as./Has) ...
              - kvas .* h0as .*exp(-(h0as./Has).^2) ...
            );
        
    Bas = C0as - Aas .* exp(-(h0as./Has).^2);
    
     concLabNew(1:indexXsmallerh) = ...
         Bas + Aas .* exp( -(xlabmesh(1:indexXsmallerh)./Has).^2 );
     fprintf('NOTE: Gaussian approx substituted for the diverging concentration h->0. ')
    
    end
     

  stateLabNew = stateLab;
    
    stateLabNew(1:indexXsmallerh) = solidStateGlobal;
    tend = t +dt; 

end

