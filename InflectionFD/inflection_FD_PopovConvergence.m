%%
close all, clear all
theta=1/2; 
WG=1;

%%
% Should follow Popov and do validation on pure WG mode?

%% Exploring dependence on dt
Z=40; dz=0.01; t0=-5; tmax=5; 
dt_vec=0.1*2.^(-(0:7));
uu=[];
n=[];
for i=1:length(dt_vec)
    dt=dt_vec(i);
    [tt,zz,u] = inflection_FD_Popov(theta,Z,t0,tmax,dz,dt,WG);    
    uu=[uu u(:,end)];
    if i>1
        n(i-1)=norm(uu(:,i)-uu(:,i-1))/norm(uu(:,i));
    end
end
figure
loglog(dt_vec(2:end),n)
xlabel('dt')
%% Exploring dependence on dz
Z=40; dt=0.05; t0=-5; tmax=5; 
dz_vec=0.1*2.^(-(7:-1:0));  % Smallest first to fix minimum mesh width
uu=[];
n=[];
for i=1:length(dz_vec)
    dz=dz_vec(i);
    [tt,zz,u] = inflection_FD_Popov(theta,Z,t0,tmax,dz,dt,WG);    
    zz_temp=zz(:,end);
    uu_temp=u(:,end);
    if i==1
        zz_fine=zz_temp;
    else
        uu_temp=interp1(zz_temp,uu_temp,zz_fine);
    end
    uu=[uu_temp uu];  % Note reversed order
end
for i=2:length(dz_vec)
    n(i-1)=norm(uu(:,i)-uu(:,i-1))/norm(uu(:,i));
end
%figure
loglog(fliplr(dz_vec(2:end)),n)
xlabel('dz')
%% Exploring dependence on t0
Z=40; dz=0.1; tmax=5; dt=0.01; 
t0_vec=[-5,-10,-15];
uu=[];
n=[];
for i=1:length(t0_vec)
    t0=t0_vec(i);
    [tt,zz,u] = inflection_FD_Popov(theta,Z,t0,tmax,dz,dt,WG);
    uu=[uu u(:,end)];
    if i>1
        n(i-1)=norm(uu(:,i)-uu(:,i-1))/norm(uu(:,i));
    end
end
%plot(abs(uu(:,2)))