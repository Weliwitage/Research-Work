function k=kernel(u)
 
%k = (1.0/sqrt(2*pi))*exp(-u*u/2);

% if ((u<1)&(u>-1)) 
%     k = (35.0/32)*(1-u*u)^3;
% else    
%   k=0.0;
% end

if ((u<1)&(u>-1)) 
    k = (3.0/4)*(1-u*u);
else    
    k=0.0;
end

% if ((u<1)&(u>-1)) 
%     k = 1/2;
%  else    
%      k=0.0;
%  end