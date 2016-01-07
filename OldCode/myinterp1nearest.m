function yi=myinterp1nearest(x,y,xi)

% for ii = 1:length(xi)
% 
%     [val,index]=min(abs(xi(ii) - x));
%     yi(ii) = y(index);
%     
% end

% 
%    [dum, j] = sort([x;xi]);
%    r(j) = 1:length(j);
%    r = r(length(x)+1:end) - (1:length(xxi));
%    r(k) = r;
%    r(xi==x(end)) = length(x)-1;
%    ind = find((r>0) & (r<length(x)));
%    ind = ind(:);
%    rind = r(ind);
%    u = (xi(ind)-xrind)./(x(rind+1)-xrind);
%    yrind = y(rind,:);
%    
%    
%    yi(ind,:)=yrind + bsxfun(@times,y(rind+1,:)-yrind,u);

   
   
   [dum, j] = sort([x;xi]);
    r(j) = 1:length(j);
    r = r(length(x)+1:end) - (1:length(xi));
    r([1:length(xi)]) = r; % Left hand index
    r(xi==x(end)) = length(x)-1; % Push end point down to other neighbour
    ind = find((r>0) & (r<length(x)));
    ind = ind(:);
    rind = r(ind);
    xrind = x(rind);
    Test = (xi(ind)-xrind)./(x(rind+1)-xrind);
    Cor = round(Test);
    yi = y(r+Cor');