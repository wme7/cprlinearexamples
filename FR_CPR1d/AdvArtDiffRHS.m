function dF = AdvArtDiffRHS(a,u,L,dg,quad,J,V,h,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Compute the Advection RHS for using NDG or FR/CPR schemes
%
%                    RHS = -f(u)_x + e(x)*u_xx
%         where f = f(x) is our Corrected/Continuous Flux function
%
%              coded by Manuel Diaz, NTU, 2014.09.16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Persson, Per-Olof, and Jaime Peraire. "Sub-cell shock capturing for
%     discontinuous Galerkin methods." AIAA paper 112 (2006): 2006. 
% [2] Kl?ckner, Andreas, Tim Warburton, and Jan S. Hesthaven. "Viscous
%     shock capturing in a time-explicit discontinuous Galerkin method."
%     Mathematical Modelling of Natural Phenomena 6.03 (2011): 57-83.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Artificial Diffusion per element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ut: Modal values
ut=V\u;

% u_hat: Truncated ut
ut(P+1,:)=0; u_hat=V*ut;

% Smooth/Resolution Indicator
s = log10(dot(u-u_hat,u-u_hat)./dot(u,u)); 

% clean 's' variable: get rid of NaN or -Inf values
s(s==-Inf)=-6; %disp(s); % considered smooth values

% Parameters
k = 1.5;	% lenght of the activation ramp
so = -4;    % $s_0$
epso = (h/P)*abs(a)/10; % $epsilon_0$

% Ranges
range1 = s<(so-k);
range2 = (so-k)<s & s<(so+k);
range3 = s>(so+k);

% Epsilon
epsilon = epso.*(0*range1 + 1/2*(1+sin(pi*(s-so)/(2*k))).*range2 + range3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 1st Derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
k1 = 0.5; % 0<=k<=1 % k=0, downwind; k=0.5, center flux, k=1; Upwind;

% u values at the boundaries of Ij (Assuming 'LGL' or 'ChebyshevMod' nodes)
switch quad
    case 'LGL'
        u_lbd = u( 1 ,:);
        u_rbd = u(end,:);
    otherwise
        u_lbd = L.lcoef*u;
        u_rbd = L.rcoef*u;
end

% Build Numerical u values across faces
u_pface = [u_lbd,0]; % + side
u_nface = [0,u_rbd]; % - side

BCs = 'Neumann';
switch BCs
    case 'Periodic' % Apply Periodic BCs
        u_nface( 1 ) = u_nface(end);    % left BD
        u_pface(end) = u_pface( 1 );    % right BD
    case 'Dirichlet' % Apply Dirichlet BCs
        u_nface( 1 ) = 0;               % left BD
        u_pface(end) = 0;               % right BD
    case 'Neumann' % Apply Neumann BCs
        u_nface( 1 ) = u_pface( 1 );	% left BD
        u_pface(end) = u_nface(end);    % right BD
end

% u communication at cell interfaces with central flux
u_comface = k1*u_nface + (1-k1)*u_pface; % common face
u_comL = u_comface(1:end-1); u_comR = u_comface(2:end);

% flux derivate
du = L.dcoef*u;

% compute auxiliary variable q 
q = (ones(P+1,1)*epsilon).*(du + dg.RR*(u_comL - u_lbd) + dg.RL*(u_comR - u_rbd))/J;

% compute fluxes boundaries of Ij
f = a*u;

switch quad
    case 'LGL'
        f_lbd = f( 1 ,:);
        f_rbd = f(end,:);
    otherwise
        f_lbd = L.lcoef*f;
        f_rbd = L.rcoef*f;
end

% LF numerical flux
alpha = abs(a);
nflux = 0.5*(a*(u_nface+u_pface)-alpha*(u_pface-u_nface));
nfluxL = nflux(1:end-1); nfluxR = nflux(2:end);

% flux derivate
df = L.dcoef*f;

% Compute the derivate: F = f + gL*(nfluxL-f_bdL) + gR*(nfluxR-f_bdR)
f_x = (df + dg.RR*(nfluxL - f_lbd) + dg.RL*(nfluxR - f_rbd))/J;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2nd Derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% switch
k2 = 0.5; % 0<=k<=1 % k=0, downwind; k=0.5, center flux, k=1; Upwind;
    
switch quad
    case 'LGL'
        q_lbd  = q( 1 ,:);
        q_rbd  = q(end,:);
    otherwise
        q_lbd  = L.lcoef*q;
        q_rbd  = L.rcoef*q;
end

% Build Numerical fluxes across faces
q_pface = [q_lbd,0]; % + side
q_nface = [0,q_rbd]; % - side

switch BCs
    case 'Periodic' % Apply Periodic BCs
        q_nface( 1 ) = q_nface(end);    % left BD
        q_pface(end) = q_pface( 1 );	% right BD
    case 'Dirichlet' % Apply Dirichlet BCs
        q_nface( 1 ) = q_pface( 1 );	% left BD
        q_pface(end) = q_nface(end);    % right BD
    case 'Neumann' % Apply Neumann BCs
        q_nface( 1 ) = 0;               % left BD
        q_pface(end) = 0;               % right BD
end

% q communication at cell interfaces with central flux
q_comface = k2*q_nface + (1-k2)*q_pface; % common face
q_comL = q_comface(1:end-1); q_comR = q_comface(2:end);

% flux derivate
dq = L.dcoef*q;

% Compute the derivate: u_xx = u_xi + gL*(nfluxL-f_bdL) + gR*(nfluxR-f_bdR)
u_xx = (dq + dg.RR*(q_comL - q_lbd) + dg.RL*(q_comR - q_rbd))/J;

% Compute RHS for the semi-discrete PDE
dF = -f_x + u_xx;