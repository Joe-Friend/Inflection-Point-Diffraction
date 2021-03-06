function [tt,xx,u] = Popov_FD_Simplified(xmax,t0,tmax,dx,dt,WG,ICorder)
% With a fixed theta = 1/2, use Popov's finite difference method to output
% a matrix of complex values from which the wave field on the domain is 
% derived.
x=(0:dx:xmax)';
M=length(x)-1;
t=(t0:dt:tmax);
N=length(t)-1;
% Uses airy0, not made by me, cited in thesis
nu = -airy0(0, WG);
%following constants for the higher order corrections
c1 = -1i * nu^2 * 32/225;
c2 = 32/45 * nu - 512/50625 * nu^4;
X = x*(-2*t0).^(1/3);
switch ICorder %Choose which initial condition order to use
    case 0
        phi0 = airy(X-nu);
        u0=(-t0)^(1/6)*exp(1i*nu*3/5*2^(-1/3)*(-t0)^(5/3))*(phi0);
    case 1
        phi0 = airy(X-nu);
        phi1 = (c1 + 1i/3 .* X.^2) .* airy(X-nu);
        u0=(-t0)^(1/6)*exp(1i*nu*3/5*2^(-1/3)*(-t0)^(5/3))*(phi0 + (-2 * t0)^(-5/3) * phi1);
    case 2
        phi0 = airy(X-nu);
        phi1 = (c1 + 1i/3 .* X.^2) .* airy(X-nu);
        phi2 = (c2 -16/45 .* X + 64/1350 * nu^2 * X.^2 -1/18 * X.^4) .* airy(X-nu) + (64/135 * nu * X + 16/45 * X.^2) .* airy(1,(X-nu));        u0=(-t0)^(1/6)*exp(1i*nu*3/5*2^(-1/3)*(-t0)^(5/3))*(phi0 + (-2 * t0)^(-5/3) * phi1 + ((-2 * t0)^(-10/3)) * phi2);
end
u = zeros(M+1,N+1);
%Set first column to the initial condition imposed
u(:,1) = u0;
%Construct the matrices required for FD scheme
LHSOffDiag = 1 / (4 * dx^2);
onesVec = ones(M-1,1);
LHSBase = spdiags([LHSOffDiag * onesVec, (1i/dt - 2 * LHSOffDiag) * onesVec , LHSOffDiag * onesVec],-1:1, M-1, M-1);
for n=1:N
    LHSFull = LHSBase + 1/2 * t(n+1) * spdiags(x(2:M),0,M-1,M-1);
    RHSFull = -conj(LHSFull);
    u(2:M,n+1)=LHSFull\(RHSFull * u(2:M,n));
end
[tt,xx]=meshgrid(t,x);