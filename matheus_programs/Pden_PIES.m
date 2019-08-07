%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TAKES Tseries AND Sseries TO CALCULATE POTENTIAL DENSITY, ASSUMING
% THAT I DON'T HAVE POTENTIAL TEMPERATURE.
%
%
%
% CORTEZI 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%close all; clear all; 
%load T_and_S_series_PIES
load Tseries;
load Sseries;
load CTD_and_Argo_p_levels

for i = 1:4
	Pdenseries(:,:,i) = sw_pden(Sseries(:,:,i),Tseries(:,:,i),dep',0);
end;

load timeaxis
clear i;
save PV_in_vars.mat 
