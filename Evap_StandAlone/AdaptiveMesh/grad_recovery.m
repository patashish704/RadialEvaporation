function dv = grad_recovery(jmax,x,v)
%
% dv = grad_recovery(jmax,x,v)
%
% this function recovers the gradient based on nodal function values.
% the derivative at a node is obtained by differentiating a quadratic
% polynomial fitting in the least squares sense the nodal function values
% at some neighboring points and evaluating the resulting polynomial at the node. 
% (x) and (v) are column vectors of size jmax.
%

%
% Copyright (C) 2010 Weizhang Huang and Robert D. Russell
% all rights reserved.
%
% This program is provided "as is", without warranty of any kind.
% Permission is hereby granted, free of charge, to use this program
% for personal, research, and education purposes. Distribution or use
% of this program for any commercial purpose is permissible
% only by direct arrangement with the copyright owner.
%

   % declare arrays

   dv = zeros(jmax,1);
   s = zeros(3,1);
   
   % main loop
   
   p = 1;
   for j=1:jmax
      % determine the range
      j1 = j-p;
      j2 = j+p;
      if j1 < 1
         j1 = 1; j2 = 2*p+1;
         if j2 > jmax
            error('jmax is too small');
         end
      end
      if j2 > jmax
         j2 = jmax; j1 = jmax-2*p;
         if j1 < 1
            error('jmax is too small');
         end
      end
      % begin the least squares fitting
      A = zeros(3,3);
      rhs = zeros(3,1);
      xc = sum(x(j1:j2))/(j2-j1+1);
      h = 1.0/max(abs(x(j1:j2)-xc));
      for k=j1:j2
         xx = (x(k)-xc)*h;
         s(1) = 1.0;
         s(2) = xx;
         s(3) = 0.5*(3*xx*xx-1.0);
         for kk=1:3, A(:,kk) = A(:,kk)+s*s(kk); end
         rhs = rhs + s*v(k);
      end
      rhs = A\rhs;
      dv(j) = (rhs(2)+3*rhs(3)*(x(j)-xc)*h)*h;
   end

% end of grad_recovery

