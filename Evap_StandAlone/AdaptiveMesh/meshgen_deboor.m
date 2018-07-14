function x = meshgen_deboor(N,y,rho,xi,M)
%MESHGEN_DEBOOR generates a new mesh of N points using de Boor's algorithm.
%
% X = MESHGEN_DEBOOR(N,Y,RHO) generates a new mesh of N points using
% de Boor's algorithm. Y is the background mesh of N points.
% RHO is the adaptation function defined on the background mesh Y. 
%
% Usage with optional arguments:
% X = MESHGEN_DEBOOR(N,Y,RHO,XI)
% X = MESHGEN_DEBOOR(N,Y,RHO,XI,M)
%
% The optional arguments are XI and M. Their default values are
% XI(i) = (i-1)/(N-1), i = 1,...,N and M = N.
%
% XI specifies the computational mesh of N points.
% 
% M is used to specify the number of points in the background mesh Y.
% M can be different from N, the number of points in the to-be-generated
% mesh. When M is not equal to N, the length of RHO must be M.
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

x = zeros(N,1);

% Initialization

if nargin < 3
   error('The first three arguments must be given.')
end
if nargin == 3
   M = N;
   xi = linspace(0,1,N)';
elseif nargin == 4
   M = N;
end

% Compute sigma

sigma = 0;
for i = 2:M
   sigma = sigma + (y(i)-y(i-1))*(rho(i)+rho(i-1));
end
sigma = sigma*0.5;

% Begin generating the new mesh

x(1) = y(1);
x(N) = y(M);
k = 2;
temp1 = 0;

for i = 2:N-1
   sigmai = xi(i)*sigma;
   for j = k:M
      temp = temp1 + 0.5*(y(j)-y(j-1))*(rho(j)+rho(j-1));
      if sigmai < temp + eps
         k = j; 
         break;
      end
      temp1 = temp;
   end 
   x(i) = y(k-1) + 2*(sigmai-temp1)/(rho(k)+rho(k-1));
end 
% end of meshgen_deboor
