% Solve the 1D beam equation with Kelvin-Voigt damping on [-1,1] with clamped BCs
% using the Galerkin method with the Chebyshev function basis introduced in
% [Shen 1995]
% Copyright (C) 2020 by Lassi Paunonen (lassi.paunonen@tuni.fi)
% Licensed under GNU GPLv3 (see LICENSE.txt).


% Basis functions are \phi_k = T_k-2*(k+2)/(k+3)*T_{k+2}+(k+1)/(k+3)*T_{k+4} for k=0..N-4, where T_k are the
% Chebyshev polynomials of the first kind

% The inner products are defined as 
% <f,g>_w = \int_{-1}^1 f(x)g(x)(1-x^2)^{-1}dx

% Compile the approximation as a system M*v''(t) + a*F*v(t)+b*F*v'(t)+c*M*v'(t)= + Bu(t)
% Single output, B=b
% The formulas for the elements M_{kj} = (<\phi_j,\phi_k>_w) 
% and A_{kj} = (<\phi_j'',(w\phi_k)''>)

% Physical parameters of the system
E = 10;
% E = 2.1e11;
I = 1;
% I = 1.167e-10;
% Damping coefficients d_KV (Kelvin-Voigt damping) and d_v (viscous
% damping)
% d_KV = .05;
d_KV = 0.01;
% d_v = 0.5;
% d_v = 0.001;
d_v = 0.4;

% Input profile functions
b1 = @(xi) 1/3*(xi+1).^2.*(1-xi).^6;
b2 = @(xi) 1/3*(xi+1).^6.*(1-xi).^2;

% Locations of the pointwise observations y(t)=[v(\xi_1,t);v(\xi_2,t)]^T
xi1 = -0.6;
xi2 = 0.3;

% Disturbance input profile function
bd1 = @(r) (r+1).^2.*(1-r).^2;
dimUd = 1; % number of disturbance signals


% Size of the higher order approximation - used for simulation, representing the original PDE
Nhi = 70; 
% Size of the lower order approximation - used for controller design
Nlo = 40;


% Definition of the reference signal y_{ref}(t) and the disturbance signal
% w_{dist}(t), and the frequencies in the signal
yref1_oneper = @(t) 0.3+0.4*((t<=1).*(t-1/2)+(t>1).*(3/2-t));
yref = @(t) [yref1_oneper(rem(t,2));zeros(size(t))];
wdist = @(t) sin(pi*t)+0.4*cos(3*pi*t);
freqs = [pi*(0:10)];
q=10;

% yref = @(t) [yref1_oneper(rem(t,2));0.1*sin(2*t)];
% wdist = @(t) sin(pi*t)+0.4*cos(3*pi*t);
% freqs = [pi*(0:10),2];

% tt = linspace(0,10);
% plot(tt,yref(tt));



% % Parameters of the controller
% Target stability margin of the closed-loop system (for LQR/LQG)
% CLstabmarg = 1;
alpha1 = 2;
alpha2 = 0.8;

Q0 = eye(2*(2*q+1));
Q1 = eye(2*(Nlo-1));
Q2 = eye(2*(Nlo-1));
R1 = eye(2);
R2 = eye(2);

% Order of the reduced order observer in the controller
ROMorder = 4;

% Parameters of the simulation
tspan = [0,16]; % time-interval of the simulation
% tspan = [0,50]; % time-interval of the simulation FOR LOW-GAIN

% Initial deflection profile (v_0) and the initial velocity (\dot{v}_0) 
% Initial state of the controller is by default zero

% v0fun = @(r) (r+1)^2.*(1-r).^3;
% v0fun = @(r) (r+1)^2.*(1-r).^2;
v0 = @(r) zeros(size(r));

v0dot = @(r) zeros(size(r));


% Parameters of the visualisation of the results
% temporal grid for the output plot
tt_output = linspace(tspan(1),tspan(2),601);

% spatial and temporal grids for the state plot
spgrid_state = linspace(-1,1,81);
tt_state = linspace(tspan(1),tspan(2),201);

% spatial and temporal grids for the animation
spgrid_anim = linspace(-1,1,130);
tt_anim = linspace(tspan(1),tspan(2),401);
anim_pause = 0.02;


% Choose whther or not to print titles of the figures
PrintFigureTitles = true;


% %% Plot the norm of the inverse \|P(is)^{-1}\| of the transfer function on iR 
% % (to study the locations of possible transmission zeros)
% 
% ss = linspace(0,20);
% Ptransfun = @(s) C*((s*eye(2*(N-1))-A)\B);
% Pinvnorms = zeros(size(ss));
% for ind = 1:length(ss)
%   tmpval=svd(Ptransfun(1i*ss(ind)));
% % tmpval=real(Ptransfun(1i*ss(ind)));
%   Pinvnorms(ind) = 1/tmpval(2);
% end
% plot(ss,Pinvnorms)


%% Construct the controller using a Galerkin approximation with size 'Nlo'

% Galerkin approximation for the controller design

Sys_Nlo = ConstrEBKVbeam(E,I,d_KV,d_v,b1,b2,xi1,xi2,bd1,Nlo);

% % Plot the eigenvalues of the beam system
% plot(eig(full(Sys_Nlo.A)),'b.','markersize',10)
% xlim([-2000,3]);

% Construct the Internal Model Based Reduced Order Controller
ContrSys = ConstrContrObsBasedROM(freqs,Sys_Nlo,alpha1,alpha2,R1,R2,Q0,Q1,Q2,ROMorder);


% % For comparison: Construct the Low-Gain Internal Model Based Controller
% % Compute the values P(iw_k) using the Galerkin approximation of the 
% % control system
% % Adjust the set of frequencies, use only k\pi for k=0,...,5
% freqs = [pi*(0:5)];
% q=5;
% 
% Pvals = cell(length(freqs),1);
% for ind = 1:length(freqs)
%   Pvals{ind} = Sys_Nlo.C*((1i*freqs(ind)*eye(size(Sys_Nlo.A))-Sys_Nlo.A)\Sys_Nlo.B);
% end
% epsgain = 0.076;
% % epsgain = [0.01,0.2]; % The algorithm can optimize CL stability margin
% [ContrSys,epsgain] = ConstrContrLGReal(freqs,Pvals,epsgain,Sys_Nlo);



%% The Closed-Loop System

% Construct the closed-loop system using a higher order Galerkin
% approximation of the control system

Sys_Nhi = ConstrEBKVbeam(E,I,d_KV,d_v,b1,b2,xi1,xi2,bd1,Nhi);

A = Sys_Nhi.A;
B = Sys_Nhi.B;
C = Sys_Nhi.C;
D = Sys_Nhi.D;
Bd = Sys_Nhi.Bd;

G1 = ContrSys.G1;
G2 = ContrSys.G2;
K = ContrSys.K;

Ae = [A,B*K;G2*C,G1+G2*D*K];
Be = [Bd, zeros(size(A,1),2);zeros(size(G1,1),dimUd),-G2];
Ce = [C,D*K];
De = -eye(2);


% Print the stability margin of the closed-loop system
Ae_evals=eig(full(Ae));
max(real(Ae_evals))

% Plot eigenvalues of the closed-loop system
figure(6)
plot(real(Ae_evals),imag(Ae_evals),'b.','markersize',15)
set(gca,'tickdir','out')
grid on
xlim([-14,0]);
% axis([-20,0,-35,35])



%% Simulation

% Define an (Nhi-1)x(Nhi+1) conversion matrix Q_coeff such that
% (c_k)_k=Q_coeff*(alpha_k)_k where c_k are the Chebyshev coefficients of a
% function and alpha_k are the coordinates of the same function in the
% basis {\phi_k}
ee = ones(Nhi-1,1);
ls = (0:(Nhi+2)).';
Dia0 = ee;
Diam2 = -2*ls(3:(end-2))./(ls(3:(end-2))+1);
Diam4 = (ls(5:end)-3)./(ls(5:end)-1);
Q_coeff = full(spdiags([Diam4,0*ee,Diam2,0*ee,Dia0],-4:0,Nhi+3,Nhi-1));


% Initial condition v0, defined as a function, 
% Chebyshev coefficients from Chebfun, converted into coordinates in \phi_k
% Express the initial function as a truncated Chebyshev expansion
v0fun = chebfun(v0,'trunc',Nhi+3);
v0dotfun = chebfun(v0dot,'trunc',Nhi+3);
% Convert the Chebyshev coefficients (coefficients of the series expansion)
% to the inner products of the function v0 with the basis functions \phi_k
v0init  = Q_coeff\chebcoeffs(v0fun);
v0dotinit  = Q_coeff\chebcoeffs(v0dotfun);

% For error checking:
% % Plot the approximation of the initial condition in the subspace V 
% v0check = chebfun(Q_coeff*v0init,'coeffs');
% plot([v0check;v0])


% Define the initial state (initial state of the controller is zero)
init = [v0init;v0dotinit];
initCL = [init;zeros(size(G1,1),1)];

% Simulate the closed-loop system
odefun = @(t,xe) Ae*xe+Be*[wdist(t);yref(t)];
sol = ode15s(odefun,tspan,initCL);



%% Plot the output 

xxe_output = deval(sol,tt_output);
yy = Ce*xxe_output;

figure(1)
hold off
clf
hold on
plot(tt_output,yref(tt_output),'--','color',0.4*[1,1,1],'linewidth',2);
plot(tt_output,yy(1,:),'color',[0, 0.4470, 0.7410],'linewidth',2);
plot(tt_output,yy(2,:),'color',[0.8500, 0.3250, 0.0980],'linewidth',2);
set(gca,'tickdir','out','PlotBoxAspectRatio',[1,0.4574,0.4574])
axis([tspan -0.2 0.55])

if PrintFigureTitles, title('The output $y(t)$ and the reference signal','Interpreter','Latex','fontsize',16), end


%% Plot the error norms \|e(t)\|
figure(2)
clf
errvals = sqrt(sum((abs(yy-yref(tt_output)).^2),1));
plot(tt_output,errvals,'color',[0, 0.4470, 0.7410],'linewidth',2);
set(gca,'tickdir','out','ytick',0:0.1:0.4,'PlotBoxAspectRatio',[1,.4006,.4006])
grid on
box off
axis([tspan -0.00 0.3])

if PrintFigureTitles, title('The tracking error norms $\|e(t)\|$','Interpreter','Latex','fontsize',16), end


%% Plot the control actions
figure(3)
uu = K*xxe_output((2*Nhi-2+1):end,:);
plot(tt_output,uu,'linewidth',2);
set(gca,'tickdir','out','PlotBoxAspectRatio',[1,.4378,.4378])
grid on
box off
ylim([-1000,1200])
if PrintFigureTitles, title('The control inputs $u_1(t)$ and $u_2(t)$','Interpreter','Latex','fontsize',16), end


%% Plot the state of the controlled beam equation
% Compute the values of the solution on the spatial grid tt_state
% by converting the coefficients (alpha_k) (given in the solution 'alphas') to
% the Chebyshev coeffients, and using them to define a Chebfun object at
% each timestep

xxvals_state = zeros(length(spgrid_state),length(tt_state));

xxe_state = deval(sol,tt_state);
Cheb_coeffs = Q_coeff*xxe_state(1:(Nhi-1),:); % Position
% Cheb_coeffs = Q_coeff*xxe(Nhi:(2*Nhi-2),:); % Velocity
for ind = 1:length(tt_state)
  cfun = chebfun(Cheb_coeffs(:,ind),'coeffs');
  xxvals_state(:,ind) = cfun(spgrid_state);
end

figure(4)
surf(tt_state,spgrid_state,xxvals_state)
set(gca,'ydir','reverse')
ylabel('$\xi$','Interpreter','Latex','fontsize',18)
xlabel('$t$','Interpreter','Latex','fontsize',18)
if PrintFigureTitles, title('The deflection of the controlled beam','Interpreter','Latex','fontsize',16), end


%% Animate the solution of the controlled beam equation
xxvals_anim = zeros(length(spgrid_anim),length(tt_anim));

xxe_anim = deval(sol,tt_anim);
Cheb_coeffs = Q_coeff*xxe_anim(1:(Nhi-1),:); % Position
% Cheb_coeffs = Q_coeff*xxe_animation(Nhi:(2*Nhi-2),:); % Velocity
for ind = 1:length(tt_anim)
  cfun = chebfun(Cheb_coeffs(:,ind),'coeffs');
  xxvals_anim(:,ind) = cfun(spgrid_anim);
end

figure(5)
% No movie recording
[~,zlims] = Animate1D(xxvals_anim,spgrid_anim,tt_anim,anim_pause,0);

% Movie recording
% [MovAnim,zlims] = Animate1D(xxvals_anim,spgrid_anim,tt_anim,anim_pause,1);

