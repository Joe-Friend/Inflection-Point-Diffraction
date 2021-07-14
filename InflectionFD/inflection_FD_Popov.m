function [tt,zz,u] = inflection_FD_Popov(theta,Z,t0,tmax,dz,dt,WG)
z=(0:dz:Z)';
M=length(z)-1
t=(t0:dt:tmax);
N=length(t)-1
if WG==1
    nu=2.338107410459763; % (Minus) first zero of Airy
elseif WG==2
    nu=4.087949444130973; % (Minus) second zero of Airy
end
u0=(-t0)^(1/6)*exp(1i*nu*3/5*2^(-1/3)*(-t0)^(5/3))*airy((-2*t0)^(1/3)*z-nu);
b = ones(M-1,1);
A = (1/(2*dz^2))*spdiags([b -2*b b], -1:1, M-1, M-1);
I=speye(M-1);
u = zeros(M+1,N+1);
u(:,1) = u0;
for n=1:N        
    Qn = t(n+1)*spdiags(z(2:M), 0, M-1, M-1);
    u(2:M,n+1)=(1i*I+theta*dt*(A+Qn))\((1i*I-(1-theta)*dt*(A+Qn))*u(2:M,n));
end
[tt,zz]=meshgrid(t,z);
end