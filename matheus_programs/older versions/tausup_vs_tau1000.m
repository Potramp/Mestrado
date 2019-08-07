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



addpath /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts
load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts/CleanUp_HydroData_ForGEM

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
p = [100:50:350];


%%%% From 100 db to 350 db
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
        

	%%%% Ploting tau(p) vs tau1000
        subplot(3,2,ip);
        str = (['plot(tau1000,tau' num2str(p(ip)) ',''k.'');']);
        eval(str);
        xlabel('tau1000');
        str=(['ylabel(''tau' num2str(p(ip)) ''');']);
        eval(str);
        str=(['title(''tau' num2str(p(ip)) ''');']);
        eval(str);

	%%%%%%%%%%%%%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
	%%%%%  CREATES VARIABLE WITH DIFF OF TAU(P) AND TAU1000;
	%%%%%  TAU(P) WILL HAVE A SECOND ROW/COLUMN FOR MONTH OF THE YEAR
	%%%%%  THIS WILL LATER BE USED TO MAKE MONTHLY BINS
	
	
	str=(['dif' num2str(p(ip)) '=NaN.*ones(length(indxf),3);']);
	eval(str);
	str=(['dif' num2str(p(ip)) '(:,1)= tau1000sc - tau' num2str(p(ip)) 'sc;']);
	eval(str);
	str=(['dif' num2str(p(ip)) '(:,2)=mm;']);
	eval(str);
	str=(['dif' num2str(p(ip)) '(:,3)=dd;']);
	eval(str);

	
end;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. POLYNOMIAL FIT -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. Plots and polyfit...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(p)
	str=(['[P S] = polyfit(tau1000sc,tau' num2str(p(i)) 'sc,3);']);
	eval(str);
	subplot(3,2,i);
	hold on;
	Y = polyval(P,tau1000sc);
	%Y = sort(Y);
	%tau1000 = sort(tau1000sc);
	plot(sort(tau1000sc),sort(Y),'m');
end;


clear temp pres mm Y P junk ip indxf indxl p ip i S str
save seasonal_correction_dif_taup_tau1000.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
