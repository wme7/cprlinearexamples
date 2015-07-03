%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Solving 1-D Advection equation with FR/CPR schemes 
%                     and Artificial Viscosity
%
%             du/dt + a*du/dx = e*u_xx,  for x \in [a,b]
%                    e = e(x): variable coefficient
%
%              coded by Manuel Diaz, NTU, 2014.01.16
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs: 
% [1] A flux reconstruction approach to high-order schemes including
%     Discontinuous Galerkin methods. H.T. Huynh, AIAA 2007.
% [2] A New Class of High-Order Energy Stable Flux Reconstruction Schemes.
%     P.E. Vincent, P. Castonguay, A. Jameson, JSC 2010.
% [3] High-Order Methods for Diffusion Equation with Energy Stable Flux
%     Reconstruction Scheme. K. Ou, P. Vincent, A Jameson, AIAA 2011-46.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Basic Scheme Implementation with SSP-RK45 intergration scheeme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
fluxfun ='linear';  % select flux function
    cfl = 0.04;     % CFL condition
      a = 20;       % advection speed 
  %tEnd = 2;        % final time
      P = 08;       % degree of accuaracy 
     nE = 080;      % number of elements
     ic = 2;        % {1} Riemann, {2} Square Jump, {3} Lifted Sine

%% PreProcess

% Build 1d mesh
xgrid = mesh1d([-0.5,1.5],nE,'LGL',P);
dx = xgrid.elementSize; J = xgrid.Jacobian; 
x = xgrid.nodeCoordinates; quad = xgrid.quadratureType;
w = xgrid.weights';	xc = xgrid.elementCenter;

% Compute gR'(xi) & gL'(xi)
RR = CorrectionPolynomial('DGRight',P); 
dg.RR = RR.eval_dP(xgrid.solutionPoints); dg.RL = -flipud(dg.RR);

% Build Lagrange k-Polynomials
l = LagrangePolynomial(xgrid.solutionPoints);
L.lcoef = double(subs(l.lagrangePolynomial,-1));
L.rcoef = double(subs(l.lagrangePolynomial,1));
L.dcoef = double(subs(l.dlagrangePolynomial,xgrid.solutionPoints));

% Vandermonde from Nodal DG tools
NDG=DGtools(xgrid.solutionPoints); V=NDG.Vandermonde2; 

% IC
switch ic
    case 1 % Riemann Problem
        u0 = IC(x,6); tEnd = 0.5;
    case 2 % Square Jump
        u0 = IC(x,8); tEnd = 1.0/abs(a);
    case 3 % Lifted Sine: 1/2 - sin(pi*x) 
        u0 = IC(x,4); tEnd = 1.5;
    otherwise 
        error('Check the list dude!')
end

% Exact solution for periodic tests
%ue = u0;

% Exact solution for Square Jump at time = t is given by,
SQJ = @(x,a,t) (heaviside(x+1/4-a*t) - heaviside(x-1/4-a*t)) + 1 ;
ue = SQJ(x,a,tEnd);

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),...
    min(min(u0))-0.1,max(max(u0))+0.1];

%% Solver Loop

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

% Using a 4rd Order 5-stage SSPRK time integration
res_u=zeros(P+1,nE); % Runge-Kutta residual storage

% Set initial time step
dt0=cfl*dx/a;

% Set initial time & load IC
t=0; u=u0; it=0; dt=dt0; 

% Build Element Error History
time=t:dt:tEnd; L1ElemErrHist=zeros(size(time,2),nE); 

while t < tEnd
    
    % Correction for final time step
    if t+dt>tEnd, dt=tEnd-t; end
    
    % iteration counter
    it=it+1;  

    % Advection Scheme with Artificial viscosity
    for RKs = 1:5
        t_local = t + rk4c(RKs)*dt;
        dF = AdvArtDiffRHS(a,u,L,dg,quad,J,V,dx,P);
        res_u = rk4a(RKs)*res_u - dt*dF;
        u = u - rk4b(RKs)*res_u;
    end

    % compute cell averages
    u_bar = w*u/2;
    
    % Update time
    t=t+dt;
      
    % Plot u
    if rem(it,10) == 0
        plot(x,u,'-+',x,ue,'--',xc,u_bar,'ro'); axis(plotrange); grid on; 
        drawnow;
    end
    
    % Element error history
    L1ElemErrHist(it,:) = dx*w*abs(u-SQJ(x,a,t));
    L1ElemErrHist(L1ElemErrHist<(1e-10)) = 1E-10;
    
end
%% Post Process

% Plot solution
figure(1)
plot(x(:),ue(:),'--k',x(:),u(:),'-+r',xc(:),u_bar(:),'ob'); axis(plotrange);
legend({'Exact solution',['FR/CPR P',num2str(P)],'Cell averages'},'Location','NorthWest','FontSize',18); legend('boxoff')
title('Scalar Advection with FR/CPR','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);

% Plot Error History
figure(2); 
imagesc(xc,time,log10(L1ElemErrHist)),h1=colorbar;set(gca,'YDir','normal');set(h1,'YDir','reverse' );
title('Convergence Error: $log_{10}(|u-u_{exact}|$)','interpreter','latex','FontSize',18);
xlabel('x','interpreter','latex','FontSize',14);
ylabel({'time'},'interpreter','latex','FontSize',14);
