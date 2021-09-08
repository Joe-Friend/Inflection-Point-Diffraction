function [tt,xx,u] = Popov_FD(theta,xmax,t0,tmax,dx,dt,WG)
x=(0:dx:xmax)';
M=length(x)-1;
t=(t0:dt:tmax);
N=length(t)-1;
if WG==1
    nu=2.338107410459763; % (Minus) first zero of Airy
elseif WG==2
    nu=4.087949444130973; % (Minus) second zero of Airy
end
u0=(-t0)^(1/6)*exp(1i*nu*3/5*2^(-1/3)*(-t0)^(5/3))*airy((-2*t0)^(1/3)*x-nu);
u = zeros(M+1,N+1);
u(:,1) = u0;
LHSOffDiag = theta / (2 * dx^2);
RHSOffDiag = (theta-1) / (2 * dx^2);
onesVec = ones(M-1,1);
LHSBase = spdiags([LHSOffDiag * onesVec, (1i/dt - 2 * LHSOffDiag) * onesVec , LHSOffDiag * onesVec],-1:1, M-1, M-1);
RHSBase = spdiags([RHSOffDiag * onesVec, (1i/dt - 2 * RHSOffDiag) * onesVec , RHSOffDiag * onesVec],-1:1, M-1, M-1);
for n=1:N
    LHSFull = LHSBase + theta * t(n+1) * spdiags(x(2:M),0,M-1,M-1);
    RHSFull = RHSBase + (theta - 1) * t(n+1) * spdiags(x(2:M),0,M-1,M-1);
    u(2:M,n+1)=LHSFull\(RHSFull * u(2:M,n));
end
[tt,xx]=meshgrid(t,x);


end