%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATES VERTICAL GRADIENT OF T ON Tseries FROM THE FILE PV_in_vars
% SO THAT A REGION OF MINIMUM GRADIENT CAN BE IDENTIFIED. THIS GRADIENT
% LAYER SHALL AFTER INDICATE A VALUE OF POTENT. VORTIC. THAT IS 
% REPRESENTATIVE OF MODE WATERS.
%
%
%
%
% -CORTEZI, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

%%%% LOADING FILE WITH Tseries FOR 200db REFERENCE LEVEL;
load ~/Dropbox/Mestrado/matheus_programs/PV_in_vars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Variables loaded. Calculating vertical gradient of T...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% PIES LOCATION
load pies_latlon;

%%%% depth
%load CTD_and_Argo_p_levels;


%%%% CALCULATING POTENTIAL TEMP.
for i = 1:4
        potTseries(:,:,i) = sw_ptmp(Sseries(:,:,i),Tseries(:,:,i),dep',0);
end;

%%%% making a matrix the same size as the temporal series but for depths(pressure)
zz = repmat(dep',[1 2429 4]);

%%%% making depth out of pressure
zz = sw_dpth(zz,pies_lat_lon(1,1));

%%%% zz is depth, so should be negative.
T_vgrad = diff(potTseries,1,1)./diff(-zz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
