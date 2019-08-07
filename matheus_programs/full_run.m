%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS WILL RUN THROUGH ALL MY USEFUL ROUTINES.
% THE IDEA IS TO TEST DIFERENT REFERENCE LEVELS FOR THE SEASONAL 
% CORRECTIONS;
%
% MUST RUN FOR T, S, CREATE VARIABLES THAT GO IN PV_PIES
% AND CALCULATE THE 16 CELSIUS ISOTHERMAL DEPTH.
% AFTER ALL THAT, COMPARE WITH PIERO (ISAS) DATA AND SELECT
% THE REFERENCE LEVEL THAT GIVES THE SMALEST DIFERENCE.
%
%
% CORTEZI, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
more('off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-.-.-.-.-.-.-.-.-. 1. Seasonal correction of tau...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tausup_vs_tau1000_v2
butt_filter;
clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-.-.-.-.-.-.-.-.-. 2. Seasonal correction of T...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
more('off');

ref_p = 100:20:360;

for w = 1:length(ref_p)
	my_p = ref_p(w);
	func_T_seasonal_model(my_p);
	func_S_seasonal_model(my_p);
end;



pause
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-.-.-.-.-.-.-.-.-. 3. Temperature time series...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
more('off');

all_p_do_all_time_series;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-.-.-.-.-.-.-.-.-. 4. 16 Celsius isothermal...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
more('off');

ref_p = 100:20:360;
for w = 1:length(ref_p)
	str = (['cd ~/Dropbox/Mestrado/matheus_programs/' num2str(ref_p(w))]);
	eval(str);
	all_p_isot16depth;
	disp('done');
	ref_p = 100:20:360;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
