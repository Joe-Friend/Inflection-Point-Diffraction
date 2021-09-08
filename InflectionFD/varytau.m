close all
clear all
figure
theta=1/2; 
WG=1;
ICorder = 1;

dtvec = [0.0005, 0.001, 0.002, 0.004, 0.008, 0.016];
Z=40; dz=0.02; tmax=5; t0=-5;
z = (0:dz:Z);
searchlightarray = cell(length(dtvec)-1, 1);
for i=1:length(dtvec)
    i
    dt = dtvec(i);
    t = (0:dt:tmax);
    [tt,zz,u] = Popov_FD_Simplified(Z,t0,tmax,dz,dt,WG,ICorder);
    [maxes, searchlightindices] = max(u);
    if i==1
        bestsearchlight = z(searchlightindices(end-(tmax/dt):end))';
        bestt = t;
    else
        tempsearchlight = z(searchlightindices(end-(tmax/dt):end))';
        searchlightarray{i-1} = tempsearchlight;
    end
end

errorarray = cell(length(dtvec)-1, 1);
tarray = cell(length(dtvec)-1,1);
for i=1:length(dtvec)-1
    tarray{i} = bestt(1:2^i:end);
    errorarray{i} = searchlightarray{i}-bestsearchlight(1:2^i:end);
end


% Recommended: run the code, then call the following lines in the command
% window with i replaced by the position required, rather than doing the
% calculations over and over again by changing the plot number here

figure;
plot(tarray{i},errorarray{i}); 
hold on;
xlabel('t');
ylabel('Difference between searchlight beams');
title('Error in searchlight beams between dt = 0.0005 and dt = [insert]');