function [xh,A,Q,cd,U,mu0,mu2,cc,gam,muc] = setup_stable(nx)

U = 2;                     % U:    Convection constant
cc = 0.2;                 % cc:   Most unstable wavenumber
cd = -1.0;                   % cd:   Dispersion

%supercrital
mu0 = -0.01;           % mu0:  Control paramter
mu2 = -.01;          % mu2:  Non-parallelity parameter

%%%% Hermite differentiation
gam = 1+sqrt(-1)*cd;			
chi = (-mu2/(2*gam))^(0.25);	% chi: decay rate of global modes
[xh, DM] = herdif(nx,2,real(chi));
D1 = DM(:,:,1);D2 = DM(:,:,2);

%%%% Intergration wieghts
w = ([diff(xh);0]+[0;diff(xh)])/2;
Q = diag(w);

%%%% Discretized Operator 
Uc = U+2*cc*real(gam)*sqrt(-1);	
mu = mu0 - real(gam)*cc^2 + mu2*xh.^2/2 ;
A = -Uc*D1+gam*D2+diag(mu);


gam = 1+sqrt(-1)*cd;
hh = sqrt(-2*mu2*gam);
Umax = U+2*cc*cd;

%%%%  Criteria for Absolute Instability and Global instability
mut=(Umax)^2/(4*abs(gam)^2);
muc=mut+abs(hh)/2*cos(angle(gam)/2);

