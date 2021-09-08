close all
clear all
figure
theta=1/2; 
WG=1;
ICorder = 1;

t0vec = [-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-20];
Z=40; dz=0.02; tmax=5; dt=0.002;
t = (t0:dt:tmax);
z = (0:dz:Z);
searchlightarray = [];
for i=1:length(t0vec)
    i
    [tt,zz,u] = Popov_FD_Simplified(Z,t0vec(i),tmax,dz,dt,WG,ICorder);
    [maxes, searchlightindices] = max(u);
    searchlightarray = [searchlightarray, z(searchlightindices(end-2500:end))'];
end
errorarray = searchlightarray - searchlightarray(:,end);

% Recommended: run the code, then call the following lines in the command
% window with T replaced by the t0 required, rather than doing the
% calculations over and over again by changing the plot number here

% figure;
% plot(t(end-2500:end),errorarray(:,T)); 
% hold on;
% xlabel('t');
% ylabel('Difference between searchlight beams');
% title('Error in searchlight beams between t_0 = -T and t_0 = -20');




