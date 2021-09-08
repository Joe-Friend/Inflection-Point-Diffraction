% Dave Hewett
% Code to plot Helmholtz solution corresponding to WG->CW inflection
% problem using an FD solution to the Popov equation.

clear all, close all
range = [-4 10 -5 10];
%range = [-10 10 -10 10];
includeRegion = @(x,y) y>-x.^3/6;
k = 80;     % wavenumber for domain plots
mu = 20;      % #pts per wavelength in domain plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=k^(-1/5)*linspace(range(1),range(2),ceil((range(2)-range(1))*mu*k^(4/5)/(2*pi))+1);
y0=k^(-3/5)*linspace(range(3),range(4),ceil((range(4)-range(3))*mu*k^(2/5)/(2*pi))+1);
[X0, Y0] = meshgrid(x0,y0);
reg = includeRegion(X0,Y0);
%X0=reg.*X0;
%Y0=reg.*Y0;
tq=k^(1/5)*reg.*X0;
zq=k^(3/5)*reg.*(Y0+X0.^3/6);
%scatter(reshape(tq,[numel(tq) 1]),reshape(zq,[numel(zq) 1]),ones(numel(tq),1))

%return
%%%%%%%%%%%%%%%%%
theta=1/2;
WG=1;
t0=min(range(1),-range(2));
tmax=range(2);
Z=tmax^3/4;
dz=0.01;
dt=0.002;
ICorder = 2;
tic
[tt,zz,u] = Popov_FD_Simplified(Z,t0,tmax,dz,dt,WG,ICorder);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[maxes, searchlightindices] = max(u);
figure
phi=exp(1i*k*(X0-.5*X0.^2.*Y0 - 7/120*X0.^5)).*(interp2(tt,zz,u,tq,zq));
%phi=exp(1i*(k*X0)).*(interp2(tt,zz,u,tq,zq));
surf(X0,Y0,real(phi),'LineStyle','none');
colorbar;
colormap jet;
view(2)
axis equal tight
hold on
%contour(X0,Y0,Y0+4*X0.^3/27,[0 0],'k','LineWidth',2)
plot3(x0,-x0.^3/6,5*ones(size(x0)),'k','LineWidth',2)
xlabel('x');
ylabel('y');
xlim([x0(1),x0(end)])
ylim([y0(1),y0(end)])

% Useful code
%print -dpng -r600 'Inflection_FirstWGMode_Helmholtz.png'
%print -dpng -r600 'Inflection_SecondWGMode_Helmholtz.png'

