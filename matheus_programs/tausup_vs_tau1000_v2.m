%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This calculates tau(p) for superficial levels and compares
% with tau1000. This is to reproduce figure A1 in Watts2001
% so a seasonal correction can be made on GEM fields.
% First, a dispersion diagram with tau(p) and tau1000 must be
% calculated.
%
% Later, a 3rd degree polynomial is fit to the curves (tau_fit).
% The residual should be tau(p) - tau_fit.
%
% 
%
%
% CORTEZI 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
more('off');


%%%%%%%% THE DATA PROVIDED BY CHRIS MEINEN IS LOADED
%%%%%%%% IT CONTAINS CTD AND ARGO PROFILES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Loading data...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% LAB
addpath /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts
load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts/CleanUp_HydroData_ForGEM

%%%% HOME
%addpath /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts
%load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts/CleanUp_HydroData_ForGEM

%%%%%% WILL BE USED LATER FOR MONTHLY BINS OF TAU;
mm = header(:,5);
dd = header(:,6);

% This gives 4 variables: header, hydro, indxl and indxf.
% The later two are the length of each cast in points and the index of the first data
% of each cast.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Calculating tau1000...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% CALCULATE tau1000
temp=hydro(:,2);
pres=hydro(:,1);
[junk,tau1000sc]=calc_tau(header,hydro,1000);


%%%% Creates a variable with the length of every profile repeating
%%%% the corresponding tau value (a single value). All profiles
%%%% are in sequence, one after the other.
%%%% I'm not sure what is the use of it.
tau1000=NaN*ones(size(temp));
for i=1:length(indxf)
  tau1000([indxf(i):indxl(i)])=tau1000sc(i);
end;


%%%%%%% CALCULATE tau(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Calculating tau(p)...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set diferent levels of p to be compared
p = [100:20:360];


%%%% From 100 db to 360 db
for ip = 1:length(p);
        temp=hydro(:,2);
        pres=hydro(:,1);

	%%%% Calculates tau
        str=(['[junk,tau' num2str(p(ip)) 'sc]=calc_tau(header,hydro,p(ip));']);
        eval(str);

	%%%% Loading variable to save computing time
        str=(['tau' num2str(p(ip)) '=NaN*ones(size(temp));']);
        eval(str);

	%%%% Creates a variable with the length of every profile repeating
	%%%% the corresponding tau value (a single value). All profiles
	%%%% are in sequence, one after the other.
	%%%% I'm not sure what is the use of it.
        for i=1:length(indxf)
          str=(['tau' num2str(p(ip)) '([indxf(i):indxl(i)])=tau' num2str(p(ip)) 'sc(i);']);
          eval(str);
        end;
        


	
end;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. POLYNOMIAL FIT -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. Plots and polyfit...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(p)
	%%%% Ploting tau(p) vs tau1000
        subplot(7,2,i);
        str = (['plot(tau1000,tau' num2str(p(i)) ',''k.'');']);
        eval(str);
        xlabel('tau1000');
        str=(['ylabel(''tau' num2str(p(i)) ''');']);
        eval(str);
        str=(['title(''tau' num2str(p(i)) ''');']);
        eval(str);

	%%%%%%%%%%%%%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
	%%%%%  CREATES VARIABLE WITH DIFF OF TAU(P) AND TAU_FIT.
	%%%%%  TAU_FIT IS A THIRD DEGREE POLYNOMIAL FIT TO THE CURVE TAU(P) vs TAU1000.
	%%%%%  TAU(P) WILL HAVE A SECOND ROW/COLUMN FOR MONTH OF THE YEAR.
	%%%%%  THIS WILL LATER BE USED TO MAKE MONTHLY BINS.
	

	%%%% Fitting polynomial curve
	str=(['[P S] = polyfit(tau1000sc,tau' num2str(p(i)) 'sc,3);']);
	eval(str);

	%%%% One subplot for each level of p used
	subplot(7,2,i);
	hold on;

	%%%% Building the polynomial
	tau_fit = polyval(P,tau1000sc);
	%Y = sort(Y);
	
	str=(['dif' num2str(p(i)) '=NaN.*ones(length(indxf),3);']);
	eval(str);


	%%%%%%%%%%%%-.-.-.-.-.-.-. CHECK WHAT IT IS THAT WORKS! -.-.-.-.-.-.-.-.-.-.-.-.-
	%%%% tau(p) - tau1000
	%str=(['dif' num2str(p(i)) '(:,1)= tau1000sc - tau' num2str(p(i)) 'sc;']);
	
	%%%% tau(p) - tau_fit
	%str=(['dif' num2str(p(i)) '(:,1)= tau' num2str(p(i)) 'sc - tau_fit;']);
	
	%%%% taufit - tau(p)
	str=(['dif' num2str(p(i)) '(:,1)= tau_fit - tau' num2str(p(i)) 'sc;']);
	%%%%%%%%%%%%-.-.-.-.-.-.-. /CHECK WHAT IT IS THAT WORKS! -.-.-.-.-.-.-.-.-.-.-.-.-
	
	eval(str);
	str=(['dif' num2str(p(i)) '(:,2)=mm;']);
	eval(str);
	str=(['dif' num2str(p(i)) '(:,3)=dd;']);
	eval(str);
	%tau1000 = sort(tau1000sc);
	
	plot(sort(tau1000sc),sort(tau_fit),'m');
end;


clear temp pres mm Y P junk ip indxf indxl p ip i S str
save ~/Dropbox/Mestrado/matheus_programs/result_variables/seasonal_correction_dif_taup_tau1000.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


