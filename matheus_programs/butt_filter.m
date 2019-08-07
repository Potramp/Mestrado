%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a first attempt to aply a Butterworth filter on
% the variables contained in the file result_variables. 
% These variables come as the result of tausup_vs_tau1000_v2.m.
% That script is used to attempt to repeat the metodology of
% Watts2001 and create a seasonal correction to the GEM 
% obtained with metodology proposed by Meinen.
% 
%
% This shall be incorporated to tausuo_vs_tau1000.m as soon as
% it's proven functional. 
%
% The goal of this part of my project is to determine a p-level
% that will represent the seasonal correction. A level below
% wich no significant addition of sazonal signal is observed.
%
% CORTEZI 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
more('off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['1. Loading all variables from tausup_vs_tau1000.m ...'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load /Users/matheusvcortezi/Dropbox/Mestrado/matheus_programs/result_variables/seasonal_correction_dif_taup_tau1000
load timeaxis;

p = [100:20:360];

month_avg_dif100 = NaN.*ones(12,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['2. Calculating monthly averages for dif(p) ...'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 1:12;
T = 12;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. Iteration -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
for i = 1:length(p);
	
	%%%% DEFINING VARIABLE TO SAVE COMPUTING TIME
	str = (['month_avg_dif' num2str(p(i)) ' = NaN.*ones(12,1);']);
	eval(str);	

	for j = 1:12
		str = (['indmon = find(dif' num2str(p(i)) '(:,2) == j);']);
		eval(str);
		str = (['month_avg_dif' num2str(p(i)) '(j) = mean(dif' num2str(p(i)) '(indmon,1));']);
		eval(str);
	end;
end;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. End of iteration -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%letra=(['rmgbyk']);
figure(1);

%%%% Plots monthly bins of averaged residuals and an annual cycle sine fit by least squares
for i = 1:6;
	subplot(3,2,i)
	hold on;
	%str = (['plot(1:12,month_avg_dif' num2str(p(i)) ',''.' letra(i) ''');']);	
	str = (['plot(1:12,month_avg_dif' num2str(p(i)) ',''.b'');']);	
	eval(str)
	title(num2str(p(i)));
	xlabel('month');
	ylabel('Residual tau(ms)');
	grid;
	str= (['A = sinfitB(1:12,month_avg_dif' num2str(p(i)) ',T);']);
	eval(str);
	newx = A(1)+A(2).*t+A(3).*sin((2*pi/T).*t)+A(4).*cos((2*pi/T).*t);
	plot(1:12,newx,'r');
	
	%%%% HIDE PV
	str = (['pv = pvar(month_avg_dif' num2str(p(i)) ',newx);']);
	
	%%%% SHOW PV
	%str = (['pv = pvar(month_avg_dif' num2str(p(i)) ',newx)']);
        eval(str);
	hold off;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['2. 3yr repmat and butter filter - 3month lowpass...']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- INTERPOLATION TO 365 DAYS -.-.-.-.-.-.-.-.-.-.-
rtau = [S_tilde S_tilde S_tilde];
days = [1:365];
mday = [31 59 90 120 151 181 212 243 273 304 334 365];
rm = [mday mday+365 mday+365*2];
rd = [days days+365 days+365*2];

%%%% loop to make the same for all depths;
for i = 1:11
        S_tilde_y(i,:) = spline(rm,rT(i,:),rd);
end;

S_tilde_y = S_tilde_y(:,366:end-365);
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- \INTERPOLATION TO 365 DAYS -.-.-.-.-.-.-.-.-.-.-

save S_tilde_y.mat S_tilde_y;
%%%% Plots the same monthly bins of residual but with a 2nd order
%%%% forwards and backwards 3-months lowpass Butterworth filter
%%%% as the curve. 
figure(2);
hold on;
for i = 1:length(p);
	str = (['rep_month_avg_dif' num2str(p(i)) ' = repmat(month_avg_dif' num2str(p(i)) ',3,1);']);
	eval(str);
	[B,A] = butter(2,0.25);
	str = (['Y(:,i) = filtfilt(B,A,rep_month_avg_dif' num2str(p(i)) ');']);
	eval(str)
	subtightplot(3,2,i,0.15,0.05,0.05)
       
	%%%% SHOW PERCENT VARIANCE ON TITLE 
	str= (['A = sinfitB(1:12,month_avg_dif' num2str(p(i)) ',T);']);
        eval(str);
        newx = A(1)+A(2).*t+A(3).*sin((2*pi/T).*t)+A(4).*cos((2*pi/T).*t);
	str=(['ylabel(''residual tau' num2str(p(i)) ''');']);
        eval(str);
       	str = (['pv = pvar(month_avg_dif' num2str(p(i)) ',newx);']);
        eval(str);
	xlabel('month');
	str=(['title('' ' num2str(p(i)) ' db - ' num2str(pv) ' percent variance'');']);
        eval(str);


	%%% Actual plot
	hold
	plot(1:13,Y(13:25,i));
	str = (['plot(1:12,month_avg_dif' num2str(p(i)) ', ''.k'' );']);
	eval(str);
	ylim([-5.1e-4 5.1e-4]);
	hold
	grid on;
end;

print -depsc2 figA1b_Watts.eps
%print -depsc2 tau_minus_taufit_figA1b_Watts.eps




%%%% Saving relevant variables, deleting uninportant ones.
clear newx dif* tau* h* month* rep_* str t dd i j pv A B T ;
save /Users/matheusvcortezi/Dropbox/Mestrado/matheus_programs/result_variables/tau_tilde.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['END OF LINE ...'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
