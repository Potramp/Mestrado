%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WE NEED T AND S TIME SERIES TO BE IN A SAME VARIABLE
% ALONG WITH OTHER USEFUL VARIABLES, LIKE LAT, LON, TIME
% AND DEPTH. THE OUTPUT HERE WILL BE USED WITH A VERSION
% OF Pden.m and then PV_*.m.
%
% - automated version to test T and PV parameters.
%
% -CORTEZI, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
more('off')

%%%%%% LOAD BOTH TIME SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Loading time series of T and S made from find_*_from_tau...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Ts;
load Ss;
%%%% pressure levels for sw_pden;
load CTD_and_Argo_p_levels

%%%% DEFINING LAST VERTICAL LEVEL;
zl = 36;


%%%%%% POTENTIAL DENSITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Computing potential density...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%%% dep HAS WRONG VALUES;
dep = dep(1:zl);

for i = 1:4
        Pdenseries(:,:,i) = sw_pden(Ss(1:zl,:,i),Ts(1:zl,:,i),dep',0);
end;

load timeaxis
clear i;
save PV_in_vars.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Potential Vorticity...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['!!!!PV_PIES!!!!'])
%PV_PIES
auto_PIES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
