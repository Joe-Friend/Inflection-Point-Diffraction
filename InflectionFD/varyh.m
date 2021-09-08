close all
clear all
figure
theta=1/2; 
WG=1;
ICorder = 1;

dzvec = [0.001, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28];
Z=40; dz=0.02; tmax=5; dt=0.002; t0=-5;
t = (t0:dt:tmax);
z = (0:dz:Z);
searchlightarray = [];
for i=1:length(dzvec)
    i
    dz = dzvec(i);
    z = (0:dz:Z);
    [tt,zz,u] = Popov_FD_Simplified(Z,t0,tmax,dz,dt,WG,ICorder);
    [maxes, searchlightindices] = max(u);
    searchlightarray = [searchlightarray, z(searchlightindices(end-2500:end))'];
end
errorarray = searchlightarray - searchlightarray(:,1);

% Recommended: run the above code, then call the following lines in the command
% window with H replaced by the h required, rather than doing the
% calculations over and over again by changing the plot number here

figure;
plot(t(end-2500:end),errorarray(:,H)); 
hold on;
xlabel('t');
ylabel('Difference between searchlight beams');
title('Error in searchlight beams between h = 0.01 and h = [insert]');


