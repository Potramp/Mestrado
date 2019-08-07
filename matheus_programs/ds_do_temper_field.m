%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a version of Chris Meinen's do_temper_field.m script, with the
% addition of a temporal dependency on tau.
% According to Watts,2001 we can derive the seasonal signal missing in
% the GEM process by folowing three steps cited in the article.
% I   - Repricate 3 times the monthly averaged residual of tau(p)* and
%       tau1000 - which is tau(p) subtracted from a 3rd deg. polynomial
%       fit to a tau(p) vs tau1000 curve.
% II  - Use a lowpass 2nd order Butterworth filter (ws = 3 months) 
%       forwards and backwards.
%
% III - Retain the curve tau_tilde for the middle year. 
%
% Finally, the deseasoned tau is: tau_ds = tau - tau_tilde.
%
%
% * I use tau(p) because the reference level is not yet defined. A number
%   of them is being tested at this stage of my work
%
% CORTEZI 2016  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% this is used after tausup_vs_tau1000_v2 and butt_filter, in that order.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['1. Loading all variables from butt_filter ...'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Load results from butt_filter
load /Users/matheusvcortezi/Dropbox/Mestrado/matheus_programs/result_variables/tau_tilde.mat
%load /Users/apple/Dropbox/Mestrado/matheus_programs/result_variables/tau_tilde.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Defining tau_tilde...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_tilde = Y(13:24,:);



%  !!!! ####  !!!! ####  !!!! ####  !!!! ####  !!!! ####  !!!! ####  !!!! ####  !!!! ####  !!!! #### !!!! #### 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. GEM in the making...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% WHAT COMES AHEAD IS A COPY OF do_temper_field.m WITH 
%%%% THE FOLLOWING CHANGES:
%%%% I  - taurange STEP IS 0.0001 INSTEAD OF 0.0010
%%%% II - tau1000sc IS CORRECTED BY TAUFIT AS IN WATTS 2001.

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. do_temper_field -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

% Script for creating the temperature field vs pressure and tau1000
%
% This script takes the input hydrographic measurements and creates 
% from them a smoothed field of temperature as a function of both 
% pressure and tau3000 (the latter is the acoustic travel time, 
% round trip, integrated between the surface and 3000 dbars).  The 
% T values are output on a regular grid of pressure and tau3000.  
%
%
%  Modified from the original GEM codes developed for the North 
%    Atlantic Current which had the following info/disclaimer.
%
%  Author:  Christopher S. Meinen, 1997
%            Graduate School of Oceanography, Univ. of Rhode Island
%
%
%  Disclaimer:  This program is provided as is without any warranty.  
%
%  This version uses 1000 dbar as the x-axis rather than 3000 dbar
%   to allow for the use of ARGO data 
% 

more off

%%%%%%%%% First load in the hydrographic data

%%%% LAB
addpath /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts
load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts/CleanUp_HydroData_ForGEM

%%%% HOME
%addpath /Users/apple/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts
%load /Users/apple/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts/CleanUp_HydroData_ForGEM


% This gives 4 variables: header, hydro, indxl and indxf.
% The later two are the length of each cast in points and the index of the first data
% of each cast.


temp=hydro(:,2);
pres=hydro(:,1);
[junk,tau1000sc]=calc_tau(header,hydro,1000);

% !!!!!! -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- CORTEZI SEASONAL CORRECTION -.-.-.-.-.-.-.-.-.-.-.-.-.-.

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- INTERPOLATION TO 365 DAYS -.-.-.-.-.-.-.-.-.-.-
rtau = [tau_tilde; tau_tilde; tau_tilde];
days = [1:365];
mday = [31 59 90 120 151 181 212 243 273 304 334 365];
rm = [mday mday+365 mday+365*2];
rd = [days days+365 days+365*2];

%%%% loop to make the same for all depths;
for i = 1:14
        tau_tilde_y(:,i) = spline(rm,rtau(:,i),rd);
end;

tau_tilde_y = tau_tilde_y(366:end-365,:);
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- \INTERPOLATION TO 365 DAYS -.-.-.-.-.-.-.-.-.-.-

%%%% THIS IS FOR USER INPUT, IN CASE IT'S NECESSARY
%disp(['p levels are ', num2str(p)]);
%p_levl = input(['Choose one p level of the above to use as seasonal correction... (1 to 6):   ']);

%%%% WE ARE USING P LEVELS FROM full_run; 
p_levl = 6;  %%%% 200dbar
%p_levl = w

pause
yy = header(:,4);
mm = header(:,5);
dd = header(:,6);


for i = 1:length(yy);
	indyear(i) = datenum(yy(i) - yy(i), mm(i), dd(i));
        if(indyear <= 365)
		%indmon = find(header(:,5) == month);
		tau1000sc(i) = tau1000sc(i) - tau_tilde_y(i,p_levl);
	else;
		tau1000sc(i) = tau1000sc(i) - tau_tilde_y(365,p_levl);
	end;
end;	


% !!!!!! -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- END OF CORTEZI SEASONAL CORRECTION -.-.-.-.-.-.-.-.-.-.-.-.-.-.


tau=NaN*ones(size(temp));
for i=1:length(indxf)
  tau([indxf(i):indxl(i)])=tau1000sc(i);
end

clear i header hydro junk %tau1000sc indxf indxl 

%%%%%%%%%%% Now build the GEM field

taurange=[1.3100:0.0001:1.3460];  % Tau grid onto which we will interpolate 
	  			  % using the spline smoothing

presrange=[0:20:1000 1050:50:5100];   % Pressure grid onto which we will 
				      % interpolate using the least-squares 
				      % splines

prange=[0:20:1000 1050:50:5100];  % Levels the hydrography is subsampled on

tempergrid=NaN*ones(length(prange),length(taurange));  % Set it up full of
  						       % NaNs to cut 
						       % computing time


% First interpolate onto regular grid of tau's

for i=1:length(prange)
 indx=find(pres==prange(i));
 tempergrid(i,:)=csaps(tau(indx),temp(indx),0.9999995,taurange);
end


%  % Next must fix the high tau end where the extrapolations get a bit 
%  % out of hand (the temps go below -2 deg C, which is not real)
%  
%  badindx=find(tempergrid<-2);
%  if ~isempty(badindx)
%    fprintf('Replacing some unrealistically low temps.\n')
%    tempergrid(badindx)=-2*ones(size(badindx));
%  end


fprintf('First interpolation is complete.\n')

% Next set up final grid, full of NaNs to cut computing time again.  

SMF_temperature=NaN*ones(length(presrange),length(taurange));

k=4;    % spline polynomial order + 1

breaks=[0 20 40:40:1000 1100:100:5000 5100]; 
   % breaks are the locations where there are knots in the temper vs
   % pressure curves.  

knots=augknt(breaks,k);

for i=1:length(taurange)
 sp=spap2(knots,k,prange,tempergrid(:,i)');
 SMF_temperature(:,i)=fnval(sp,presrange)';
end

%  % Again must fix the high tau end where the extrapolations get a bit 
%  % out of hand (the temps go below -2 deg C, which is not real)
%  
%  badindx=find(SMF_temperature<-2);
%  if ~isempty(badindx)
%    fprintf('Replacing some unrealistically low temps again.\n')
%    SMF_temperature(badindx)=-2*ones(size(badindx));
%  end
%  clear badindx

fprintf('Second interpolation is complete.\n')

%clear breaks k knots sp


%************************************************************************
% Error map stuff
%
% The next section of the script calculates error information for the 
% resulting smoothed field of temperature.
%

errtaubin=[1.3135:0.003:1.3435];  % Set up tau bins for calculating errors 
errpresbin=[50:100:5100];         % Set up pressure bins for calculating errors

% Set up error result arrays.  One set for absolute errors and
% one set for errors as a percentage of signal at a given depth.

biasgridpct=NaN*ones(length(errpresbin),length(errtaubin));
biasgridabs=NaN*ones(length(errpresbin),length(errtaubin));
scattergridpct=NaN*ones(length(errpresbin),length(errtaubin));
scattergridabs=NaN*ones(length(errpresbin),length(errtaubin));
rmsgridpct=NaN*ones(length(errpresbin),length(errtaubin));
rmsgridabs=NaN*ones(length(errpresbin),length(errtaubin));

% Determine the temperature values predicted by the smoothed field at 
% given the pressure and tau values of the original hydrography.    

predictedT=interp2(taurange,presrange,SMF_temperature,tau,pres);

% Finally, determine the differences between the real measured 
% temperature values and the temperature values predicted by the 
% smoothed field.  

for i=1:length(errpresbin)
 indx=find(pres>=(errpresbin(i)-50) & pres<=(errpresbin(i)+50));
 gindx=find(presrange>=(errpresbin(i)-50) & presrange<=(errpresbin(i)+50));
  temprange=max(mean(SMF_temperature(gindx,:))) - ...
            min(mean(SMF_temperature(gindx,:)));
 for j=1:length(errtaubin)
  indx2=find(tau(indx)>=(errtaubin(j)-0.0015) & ... 
                                      tau(indx)<=(errtaubin(j)+0.0015));
  gindx2=find(taurange>=(errtaubin(j)-0.0015) & ... 
		                taurange<=(errtaubin(j)+0.0015));
  if length(indx2) >= 12
  %  if ~isempty(indx2)
    biasgridpct(i,j)=100*(mean(temp(indx(indx2))) - ... 
                   mean(mean(SMF_temperature(gindx,gindx2))))/temprange;
    scattergridpct(i,j)=100*(mean(std(temp(indx(indx2)) - ...
                     predictedT(indx(indx2)))))/temprange;
    rmsgridpct(i,j)=100*(mean(rms(temp(indx(indx2)) - ...
                     predictedT(indx(indx2)))))/temprange;
    biasgridabs(i,j)=mean(temp(indx(indx2))) - ... 
                                mean(mean(SMF_temperature(gindx,gindx2)));
    scattergridabs(i,j)=mean(std(temp(indx(indx2)) - predictedT(indx(indx2))));
  rmsgridabs(i,j)=mean(rms(temp(indx(indx2)) - predictedT(indx(indx2))));

  end
 end
end

clear indx gindx indx2 gindx2 i j temprange

fprintf('Error calculations are complete.\n')

%************************************************************************

more on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


