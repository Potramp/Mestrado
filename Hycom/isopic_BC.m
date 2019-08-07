%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TO DETERMINE THE CORRECT DEPTH OF THE LOWER LIMIT OF THE BC,
% WE CHOSE THE ISOPYCNAL OF 32.2 SYGMA1000. 
% HERE THE AVERAGE ISOPYCNALS ARE CALCULATED.
%
%
%
% - CORTEZI, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
more('off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Loading data...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ~/Dropbox/Mestrado/Hycom/hycom_34S.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Potential desity...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% AVERAGING SALINITY AND TEMPERATURE FROM HYCOM DATA
s = nanmean(sal,3);
t = nanmean(tem,3);

%%%% POTENTIAL DENSITY
pp = repmat(dep',[1 101]);
pd = sw_pden(s,t,pp,1000);

%%%% PLOT
contourf(lon,-dep,pd-1000,30)                 
hold on
contour(lon,-dep,pd-1000,[32.2 32.2],'Color',[0 0 0]+0.5,'linewidth',2);
colorbar
title('Isopycnal of sigma theta(1000) equal 32.2');
xlabel('Longitude');
ylabel('Depth');
print -depsc isopicCB.eps

