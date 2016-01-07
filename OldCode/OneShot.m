function [concLabNew,stateLabNew,tend,hverst,Dtsaved] = OneShot(hIC,TIC, concLabIC, concLabbulkIC,stateLabIC, stateLabbulkIC,tIC)

global K1 K2 Z xmesh plotflag
global dt dtmax dhmax convNewton convthreshld maxIterNewton maxSteps deltaTsave
global xlabmesh 
global solidStateGlobal concSolidEnforceGlobal
global concLab concLabpp stateLab tglobal concLabbulk stateLabbulk 

%Niall's global
%global DTOutSave


concLab = concLabIC;
concLabbulk = concLabbulkIC;
stateLab = stateLabIC;
stateLabbulk = stateLabbulkIC;

solidStateIC = getInterfaceStateFromstateLab(hIC); % at the interface

hICinOneShot = hIC
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
hbreak = 0.5*hIC;  % value of h where we end the computation
%hbreak = 0.99*hIC;  % value of h where we end the computation
% This is a globally defined variable that contains the (interpolated)
% concentration in the solid at position x=h (enforced during the melting
% process but NOT during resolidification)
concSolidEnforceGlobal = CIC(K1+1);
% Global variable for the solid state AT the interface
solidStateGlobal = solidStateIC;
% add time as a globally available quantity (for explicitly time dependent
% conditions)
tglobal =t;

uold=uIC;

solidflag = false; % indicate that we are resolidifying (hdot <0)


% for output purposes
% hverst = zeros(500,2); % preaccolacte memory
% Dtsaved = zeros(500,2);
tsaveout = t;
savecounter = 1;

for ii = 1:maxSteps % MAIN LOOP
      %  ii,t
     % variable timestep version similar to Scotts code
        % first we need to take an accurate step
    accurate = false ;      
    

    while (~accurate)
        % take two dt/2-steps, then a dt-step
        success = true ;
         uoldsize = size(uold);
        if (success) ; [guess2, success] = IFunc(uold, dt/2) ; end ;
        if (success) ; [guess2, success] = IFunc(guess2, dt/2) ; end ;
        if (success) ; [guess1, success] = IFunc(uold, dt, guess2) ; end ;

        % if successful, then compare for error
        if (success) ;
           
            % CHECK NORM USED FOR ERROR ESTIMATES - IT"S A LOCAL ONE!!!!
          %  err = max( (guess2 - guess1)./(guess2+convthreshld));
            err = max( (guess2 - guess1));
            accurate = (err < convthreshld); %*dt/tfinal) 

            % Initical Temperature profile has jumps... maybe reduce the
            % accuracy requirements
            if (ii <50) ; accurate = (err < 2e-1) ; end;
            if (ii <15) ; accurate = (err < 2) ; end;
        end ;

        if (~success)
                
                fprintf('Newton not converging\n')
        end
        % cut timestep if unsuccessful or inaccurate
        if (~accurate)
            dt = dt / 2 ;
            if (dt < 1e-17) ;  
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
        concSolidEnforceGlobal = getInterfaceCSfromCSlab(h);
    
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
        
        % set all grid points in the lab frame < h to the concentration at x=h 
        % or interpolate after moving past the point !
        concLab(1:indexXsmallerh) = solidconc;
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
    plot(xmesh(1:K1+ceil(2*K2/4))*h,guess2(1:K1+ceil(2*K2/4)),'r'); % Temperature
    title('Temp');
           axis([0 10 1 5])
    %hold on;
    subplot(4,1,2)
    plot(xmesh(1:K1+1)*h,guess2(Z+1:Z+K1+1),'b'); % concentrtion in liquid
    title('concentration (liquid)');
    subplot(4,1,3)
    plot(xlabmesh,concLab,'b'); % concentration in the solid
    title('concentration (solid)')
    drawnow;    
    end

    
    if (t+dt>tsaveout+ deltaTsave)
        hverst(savecounter,1) = t+dt;
        hverst(savecounter,2) = guess2(end);
%         Dtsaved(savecounter,1) = t+dt;
%         DTOutSave= DTOutSave;
%         Dtsaved(savecounter,2) = max(DTOutSave);
        tsaveout = t+dt;
        savecounter=savecounter +1;
       
%         if (plotflag > 0)
%         figure(5) 
%         plot(hverst(:,1), hverst(:,2))
%             drawnow; 
       % end
       
    end
    
    
    
    
    
    
    if guess2(end) <hbreak  
        break
    end
    
    uold = guess2;
    t = t+dt;
    tglobal = t;

   
    
end

% Return values

concLabNew = concLab;
% Special treatment for the points smaller than hbreak:
     indexXsmallerh = find(xlabmesh<h,1,'last');
     indexXlargerh  = find(xlabmesh>h,1,'first');
     
    concLabNew(1:indexXsmallerh) = ...
         interp1([h,xlabmesh(indexXlargerh+1)],[solidconc,concLab(indexXlargerh+1)],xlabmesh(1:indexXsmallerh),'linear','extrap');
    stateLabNew = stateLab;
    
    stateLabNew(1:indexXsmallerh) = solidStateGlobal;
    tend = t -dt; 

end

