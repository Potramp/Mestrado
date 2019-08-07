% Script for creating the salinity field vs pressure and tau1000
%
% This script takes the input hydrographic measurements and creates 
% from them a smoothed field of salinity as a function of both 
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

addpath /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts
load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts/CleanUp_HydroData_ForGEM

% This gives 4 variables: header, hydro, indxl and indxf.
% The later two are the length of each cast in points and the index of the first data
% of each cast.


sal=hydro(:,3);
pres=hydro(:,1);
[junk,tau1000sc]=calc_tau(header,hydro,1000);

tau=NaN*ones(size(sal));
for i=1:length(indxf)
  tau([indxf(i):indxl(i)])=tau1000sc(i);
end

clear i header hydro indxf indxl tau1000sc junk

%%%%%%%%%%% Now build the GEM field

taurange=[1.3100:0.0001:1.3460];  % Tau grid onto which we will interpolate 
	  			  % using the spline smoothing

presrange=[0:20:1000 1050:50:5100];   % Pressure grid onto which we will 
				      % interpolate using the least-squares 
				      % splines

prange=[0:20:1000 1050:50:5100];  % Levels the hydrography is subsampled on

salgrid=NaN*ones(length(prange),length(taurange));  % Set it up full of
  						       % NaNs to cut 
						       % computing time


% First interpolate onto regular grid of tau's

for i=1:length(prange)
 indx=find(pres==prange(i));
 salgrid(i,:)=csaps(tau(indx),sal(indx),0.9999995,taurange);
end


%  % Next must fix the high tau end where the extrapolations get a bit 
%  % out of hand (the sals go below -2 deg C, which is not real)
%  
%  badindx=find(salgrid<-2);
%  if ~isempty(badindx)
%    fprintf('Replacing some unrealistically low sals.\n')
%    salgrid(badindx)=-2*ones(size(badindx));
%  end


fprintf('First interpolation is complete.\n')

% Next set up final grid, full of NaNs to cut computing time again.  

SMF_salinity=NaN*ones(length(presrange),length(taurange));

k=4;    % spline polynomial order + 1

breaks=[0 20 40:40:1000 1100:100:5000 5100]; 
   % breaks are the locations where there are knots in the sal vs
   % pressure curves.  

knots=augknt(breaks,k);

for i=1:length(taurange)
 sp=spap2(knots,k,prange,salgrid(:,i)');
 SMF_salinity(:,i)=fnval(sp,presrange)';
end

%  % Again must fix the high tau end where the extrapolations get a bit 
%  % out of hand (the sals go below -2 deg C, which is not real)
%  
%  badindx=find(SMF_salinity<-2);
%  if ~isempty(badindx)
%    fprintf('Replacing some unrealistically low sals again.\n')
%    SMF_salinity(badindx)=-2*ones(size(badindx));
%  end
%  clear badindx

fprintf('Second interpolation is complete.\n')

%clear breaks k knots sp


%************************************************************************
% Error map stuff
%
% The next section of the script calculates error information for the 
% resulting smoothed field of salinity.
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

% Determine the salinity values predicted by the smoothed field at 
% given the pressure and tau values of the original hydrography.    

predictedS=interp2(taurange,presrange,SMF_salinity,tau,pres);

% Finally, determine the differences between the real measured 
% salinity values and the salinity values predicted by the 
% smoothed field.  

for i=1:length(errpresbin)
 indx=find(pres>=(errpresbin(i)-50) & pres<=(errpresbin(i)+50));
 gindx=find(presrange>=(errpresbin(i)-50) & presrange<=(errpresbin(i)+50));
  salrange=max(mean(SMF_salinity(gindx,:))) - ...
            min(mean(SMF_salinity(gindx,:)));
 for j=1:length(errtaubin)
  indx2=find(tau(indx)>=(errtaubin(j)-0.0015) & ... 
                                      tau(indx)<=(errtaubin(j)+0.0015));
  gindx2=find(taurange>=(errtaubin(j)-0.0015) & ... 
		                taurange<=(errtaubin(j)+0.0015));
  if length(indx2) >= 12
  %  if ~isempty(indx2)
    biasgridpct(i,j)=100*(mean(sal(indx(indx2))) - ... 
                   mean(mean(SMF_salinity(gindx,gindx2))))/salrange;
    scattergridpct(i,j)=100*(mean(std(sal(indx(indx2)) - ...
                     predictedS(indx(indx2)))))/salrange;
    rmsgridpct(i,j)=100*(mean(rms(sal(indx(indx2)) - ...
                     predictedS(indx(indx2)))))/salrange;
    biasgridabs(i,j)=mean(sal(indx(indx2))) - ... 
                                mean(mean(SMF_salinity(gindx,gindx2)));
    scattergridabs(i,j)=mean(std(sal(indx(indx2)) - predictedS(indx(indx2))));
  rmsgridabs(i,j)=mean(rms(sal(indx(indx2)) - predictedS(indx(indx2))));

  end
 end
end

clear indx gindx indx2 gindx2 i j salrange

fprintf('Error calculations are complete.\n')

%************************************************************************

more on
