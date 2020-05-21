function Sys = ConstrEBKVbeam(E,I,d_KV,d_v,b1,b2,xi1,xi2,bd1,N)
% Construct a spectral Galerkin approximation with the Chebyshev function 
% basis functions in [Shen 1995] for the 1D beam equation with Kelvin-Voigt 
% damping on [-1,1] with clamped boundary conditions at both ends.
% The system has two inputs with profile functions 'b1' and 'b2', 
% two measured outputs, which are the deflections of the beam at 'xi1' and 
% 'xi2', and a single disturbance input with a profile function 'bd1'. The
% size of the resulting approximate system is (N-1)x(N-1)
%
% Copyright (C) 2020 by Lassi Paunonen (lassi.paunonen@tuni.fi)
% Licensed under GNU GPLv3 (see LICENSE.txt).


ee = ones(N-1,1);
kk = (0:(N-2))';

Dia0 = pi/2*([2;ones(N-2,1)]+(4*(kk+2).^2+(kk+1).^2)./(kk+3).^2);
Dia2 = -pi*((kk+2)./(kk+3)+(kk+1).*(kk+4)./(kk+3)./(kk+5));
Dia4 = pi/2*(kk+1)./(kk+3);

% The mass matrix for the second order system
M = spdiags([Dia4,0*ee,Dia2,0*ee,Dia0,0*ee,[0;0;Dia2(1:(end-2))],0*ee,[zeros(4,1);Dia4(1:(end-4))]],-4:4,N-1,N-1);

% Construct the matrix A
F = zeros(N-1,N-1);
for ind = 1:(N-1)
  l = ind-1;
  F(ind,ind) = 8*pi*(l+1).^2*(l+2)*(l+4);
  if ind<N-1
    ks = (l+2):2:(N-2);
    F(ind,(ind+2):2:end) = 8*pi*(l+1)*(l+2)./(ks+3).*(l*(l+4)+3*(ks+2).^2);
  end
end

A0 = M\F;
A = [zeros(N-1),eye(N-1,N-1);-E*I*A0,-d_KV*I*A0-d_v*eye(N-1)];

% Construct the input matrix B

% The corresponding input matrix B is a vector of the inner products
% <b_j,\phi_l> for l=0..N-2. These coefficients are computed from the 
% Chebyshev coefficients of the functions b_j(.), using ||T_0||_w^2= pi and 
% ||T_l||_w^2= pi/2 for l>0

b1fun = chebfun(b1,'trunc',N+3);
b2fun = chebfun(b2,'trunc',N+3);
% plot([b1fun,b2fun])

% Compute the Chebyshev inner products <b,\phi_l>_w of b
ee = ones(N-1,1);
kk = (0:(N-2)).';
Dia0 = ee;
Dia2 = -2*(kk+2)./(kk+3);
Dia4 = (kk+1)./(kk+3);

QB=full(spdiags([Dia0,0*ee,Dia2,0*ee,Dia4],0:4,N-1,N+3));

B0 = M\(QB*diag([pi;pi/2*ones(N+2,1)])*[chebcoeffs(b1fun),chebcoeffs(b2fun)]);
B = [zeros(N-1,2);B0];

% Construct the output matrix C
C_polyT= chebpoly(0:(N+2));
C0 = [C_polyT(xi1);C_polyT(xi2)]*QB';
C = [C0,zeros(2,N-1)];


% The system has zero feedthrough
D = 0;

% Construct the disturbance input matrix Bd
bd1fun = chebfun(bd1,'trunc',N+3);
dimUd = 1; % number of disturbance signals
Bd0 = M\(QB*diag([pi;pi/2*ones(N+2,1)])*chebcoeffs(bd1fun));
Bd = [zeros(N-1,dimUd);Bd0];


Sys.A = A;
Sys.B = B;
Sys.C = C;
Sys.D = D;
Sys.Bd = Bd;

