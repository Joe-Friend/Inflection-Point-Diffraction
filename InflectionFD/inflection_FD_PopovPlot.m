close all
clear all
figure
theta=1/2; 
WG=1;
ICorder = 2;

% Z=12; dz=0.03; t0=-7; tmax=3; dt=0.05; % Popov/Pschenchik parameters (I set
% Z=12 as they said they used h=0.03 and a maximum of M=400 spatial elements)

% Z=250; dz=0.03; t0=-10; tmax=15; dt=0.001; 

Z=175; dz=0.01; t0=-4; tmax=10; dt=0.002;
%Z=40; dz=0.05; t0=-15; tmax=5; dt=0.001;
%Z=250; t0=-10; tmax=10; dz=0.01; dt=0.002;
tic
[tt,zz,u] = Popov_FD_Simplified(Z,t0,tmax,dz,dt,WG,ICorder);
toc
[maxes, searchlightindices] = max(u);
Np=301;
tp = linspace(t0,tmax,Np);
zp = linspace(0,Z,Np);
[ttp,zzp]=meshgrid(tp,zp);
up = interp2(tt,zz,u,ttp,zzp);
%surf(tt,xx,real(u),'LineStyle','none')
surf(ttp,zzp,abs(up),'LineStyle','none', 'DisplayName', '\psi(x,t)')
hold on
tplus=linspace(0,tmax,Np);
plot3(tplus,tplus.^3/6,ones(size(tplus)),'b','LineWidth',1, 'DisplayName', 'x=t^3/6')
t = (t0:dt:tmax);
z = (0:dz:Z);
plot3(t,searchlightindices*dz,ones(size(t)),'g','LineWidth',1, 'DisplayName', 'Searchlight')
view(2)
xlim([t0, tmax]);
colormap hot
colorbar;
%caxis([-0.2,1])
%axis square
xlabel('t'), ylabel('x')
legend('Location', 'northwest');
filename=['Z' num2str(Z) '_dz' num2str(dz) '_t0' num2str(t0) '_tmax' num2str(tmax) '_dt' num2str(dt)];
if WG==1
    % title('First WG mode')
    % print -dpng -r300 'Inflection_FirstWGMode.png'
elseif WG==2
    % title('Second WG mode')
    % print -dpng -r300 'Inflection_SecondWGMode.png'
end