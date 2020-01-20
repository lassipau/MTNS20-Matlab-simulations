function ContrSys = ConstrContrObsBasedROM(freqs,Sys,alpha1,alpha2,R1,R2,Q0,Q1,Q2,ROMorder)
% ContrSys = ConstrContrObsBasedReal(freqs,Pvals,Sys)
%
% Construct an observer-based robust controller for systems with the same number of 
% inputs and outputs. The frequencies are assumed to be conjugate pairs, and the internal 
% model is in real form
% freqs = Frequencies to be included in the controller, only real nonnegative
% frequencies, if zero frequency is included, it's the first element in the
% vector. The control system is assumed to be real (i.e.,
% P(conj(s))=conj(P(s))), and P(iw_k) are invertible at the frequencies of
% the reference and disturbance signals.
% Sys = The Galerkin approximation (A_N,B_N,C_N,D) of the control system 
% for the controller design
% IMstabmarg = intended stability margin of the closed-loop system
% ROMorder = order of the reduced-order observer in the controller
%
% ContrSys = Controller parameters (ContrSys.G1,ContrSys.G2,ContrSys.K)
%
% In this version, the parameters in the LQR/LQG design Q are chosen to be
% identity matrices.
%
% Copyright (C) 2020 by Lassi Paunonen (lassi.paunonen@tuni.fi)
% Licensed under GNU GPLv3 (see LICENSE.txt).

A = Sys.A;
B = Sys.B;
C = Sys.C;
D = Sys.D;

dimX = size(A,1);
dimY = size(C,1);
dimU = size(B,2);

if dimY ~= dimU
  error('The has an unequal number of inputs and outputs, the observer-based controller design cannot be completed.')
end

q = length(freqs);


if freqs(1)==0, dimZ = dimY*(2*q-1); else dimZ = dimY*2*q; end

% B1 = zeros(dimZ,dimY);
% if freqs(1)==0
%   offset = 1;
%   B1(1:dimY,:) = PKvals;
% end
% 
% 
% PKappr = @(s) (C+D*K21)*((s*eye(dimX)-(A+B*K21))\B)+D;
% for ind = 1:q
%   if cond(PKappr(1i*freqs(ind)))>1e6
%     warning(['The matrix P_K(iw_k) for k=' num2str(ind) ' is nearly singular!'])
%   end
% end


% Construct the internal model
[G1,G2] = ConstrIMReal(freqs,dimY);

dimZ = size(G1,1);

As = [G1,G2*C;zeros(dimX,dimZ),A];
Bs = [G2*D;B];

Qs = blkdiag(Q0,Q2);

SigmaN = are((A+alpha1*eye(dimX))',C'*(R1\C),Q1*Q1');
PiN = are(As+alpha2*eye(dimZ+dimX),Bs*(R2\Bs'),Qs'*Qs);

L = -SigmaN*(C'/R1);
K = - (R2\Bs')*PiN;
K1N = K(:,1:dimZ);
K2N = K(:,(dimZ+1):end);

rsys = balred(ss(A+L*C,[B+L*D,L],K2N,zeros(dimU,dimU+dimY)),ROMorder);
ALr = rsys.A;
Br_full = rsys.B;
BLr = Br_full(:,1:dimU);
Lr = Br_full(:,(dimU+1):end);
K2r = rsys.C;

% % Debug: Omit model reduction
% ALr = A+L*C;
% BLr = B+L*D;
% Lr = L;
% K2r = K2N;
% ROMorder = dimX;

ContrSys.G1 = [G1,zeros(dimZ,ROMorder);BLr*K1N,ALr+BLr*K2r];
ContrSys.G2 = [G2;-Lr];
ContrSys.K = [K1N, K2r];



% % Temp
% % Check validity of the controller!
% 
% Gf1 = ContrSys.G1;
% Gf2 = ContrSys.G2;
% Kf = ContrSys.K;
% Ae = [A B*Kf;Gf2*C Gf1+Gf2*D*Kf];
% Qe = [-eye(dimX) zeros(dimX,size(Gf1,2));H zeros(size(G1,1),size(Gf1,2)); -eye(dimX) zeros(dimX,size(G1,1)) eye(dimX)];


