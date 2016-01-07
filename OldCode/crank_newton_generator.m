function ff = crank_newton_generator(residual, jacobian,convNewton,maxIter)
global plotflag jacob
 %plotflag =2;
% iter_func = crank_newton_generator(residual, jacobian)
% -------------------------------------------------------------------------
% in the context of a generic differential equation solver, this function 
% takes as arguments the function handles @residual and @jacobian, where
%
%    [resid] = residual(uold, uguess, dt),
%
% returns a vector, and
%
%    [jacob] = jacobian(uold, uguess, dt),
%
% return a matrix.  It then returns a handle to an *iterator function*, 
% @iterfunc, with an interface:
% 
%    [endvals, success] = iter_func(startvals, dt).
% 
% The resulting function @iterfunc takes as arguments a vector $startvals 
% of current function values, and the timestep $dt that you want to advance.
% It then attempts to find the value in the future, and indicates whether
% it succeeded or not.
%


% return the prepared function
ff = @nested_newton ;


function [endvals success] = nested_newton(startvals,dt,guess)

  % our initial guess for endvals will simply be startvals
  if (nargin < 3) ; guess = startvals ; end ;

  % Calculate the residual using the user-supplied function.
  resid = feval(residual, startvals, guess, dt) ;
  
  % Now we improve the guess iteratively using Newton's method.  Since
  %     (guess - actual) \approx JACOB * resid,
  % we calculate the jacobian using the user-supplied variation,
  % invert it, multiply by the residual and subtract from the guess
  % to obtain an improved guess.  We do this repeatedly until the
  % residual of the guess converges to zero.  We want a very low 
  % residual error to ensure stability.
  newtons = 0 ;
  success = false ;
  abort   = false ;
  while (success == false)
    
    % try a newton step and see if we converged
    lastwarn('') ;
    jacob = feval(jacobian, startvals, guess, dt) ;  
    
    %if(condest(jacob)>1e6) 
    %    condest(jacob)
    %end
    %full(jacob(:,:))
    %HeatMap((abs(flipud(jacob))))

    guess = guess - jacob\resid ; 
    
    %error
    
    resid = feval(residual, startvals, guess, dt) ;
    if (plotflag>1)
    subplot(4,1,4)
    plot(abs(resid)); title('abs(resid)'); drawnow;
    end
    rnorm = max(abs(resid)) ;
    success = (rnorm < convNewton) ;
    newtons = newtons + 1 ;
    
   %     tempcheck(newtons) = rnorm;
    
    
    % if the step didn't work or we are not converging, abort
    if (success==false && (~isempty(lastwarn) || newtons > maxIter))
      success = false ;
      endvals = guess ;
      return ;
    end
  end
  endvals = guess ;
  
%  %   subplot(4,1,4)
   
%   
%   if newtons > 20
%      error 
%   end
end

end
