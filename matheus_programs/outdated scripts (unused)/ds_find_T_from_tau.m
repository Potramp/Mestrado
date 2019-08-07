% THIS CODE PLOTS A TIME SERIES FOR THE YEAR 2009 FROM THE DATA IN SAM_PIES_data.txt
% USING THE GEM FIELD DONE BY ds_do_temper_field.m

% CREATES GEM
ds_do_temper_field

%LOADS PIES DATA
load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/PIES_data/SAM_PIES_data.txt
%load /Users/apple/Dropbox/Mestrado/chris_programs/PIES_data/SAM_PIES_data.txt

% TAU DATA FOR PIES-D
tau_pies = SAM_PIES_data(:,[4 6 8 10]);

% DATE VECTOR (yy,mm,dd)
date_pies = SAM_PIES_data(:,1:3);
numdate_pies = datenum(date_pies);
mm = date_pies(:,2);
% T_tilde FOR SEASONAL CORRECTION
load T_tilde

%%%% THIS SECTION IS NOT BEING USED ANYMORE. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- FOR 2009 ONLY TIME SERIES -.-.-.-.-.-.-.-.-.-.-.-
%% SELECTING ONLY THE FIRST YEAR
%ind2009 = find(date_pies(:,1) == 2009);
%tau2009 = tau_pies(ind2009);

%% LOADING VARIABLE TO SAVE COMPUTING TIME
%Tseries = NaN.*ones(133,length(tau2009));

% FINDS THE SMF_Temperature PROFILE THAT HAS THE TAU VALUE CLOSEST TO THAT OF THE
% PIES DATA
%for i = 1:length(tau2009);
%	similar = find(taurange<=tau_pies(i)+0.001 & taurange >= tau_pies(i)-0.001);
%	tam = size(similar);
%	if(tam(2)>0)
%		closest = (tau2009(i)-taurange(similar));
%		closest = abs(closest);
%		closest = (closest == min(closest));
%		indclose = find(closest == 1);
%		Tseries(:,i) = SMF_temperature(:,similar(indclose));
%	else
%		Tseries(:,i) = NaN.*ones(1,133);
%	end;
%end;
%timeaxis = numdate_pies(1:289)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- FOR ALL TIME SERIES -.-.-.-.-.-.-.-.-.-.-.-

Tsize = size(tau_pies);

% LOADING VARIABLE TO SAVE COMPUTING TIME
Tseries = 0.*ones(133,Tsize(1));


%%% PREPARING FOR SEASONAL CORRECTION
%%% Making T_tilde the size as SMF_temperature filling with NaNs.
junk = size(SMF_temperature);
junk = 0.*ones(junk(1),12);
junk2 = size(T_tilde);
junk(1:junk2,:) = T_tilde;
T_tilde = junk;

clear junk*; 

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. RECONSTRUCTING TIME SERIES -.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reconstructing time series of T  ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% FINDS THE SMF_Temperature PROFILE THAT HAS THE TAU VALUE CLOSEST TO THAT OF THE
% PIES DATA
for i = 1:length(tau_pies);
       similar = find(taurange<=tau_pies(i)+0.001 & taurange >= tau_pies(i)-0.001);
       tam = size(similar);
       if(tam(2)>0)
               closest = (tau_pies(i)-taurange(similar));
               closest = abs(closest);
               closest = (closest == min(closest));
               indclose = find(closest == 1);
		%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. SEASONAL CORRECTION -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
		%%% I need to make the seasonal correction for T as in Watts 2001.
               Tseries(:,i) = SMF_temperature(:,similar(indclose)) + T_tilde(:,mm(i));
       else
               Tseries(:,i) = NaN.*ones(1,133);
       end;
end;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. RECONSTRUCTING TIME SERIES -.-.-.-.-.-.-.-.-.-.-.-.-.

timeaxis = numdate_pies;

%%%%-.-.-.-.-.-.-.-.-.-.-.-.-.  Calculating Tg -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%% For T seasonal cycle as in Watts2001, we apply the same method did in ds_do_Temper_field
%%%% for tau, but this time for T. Tg is calculated from the tau data derived from CTD.
%%%% T_prime  will be the annual march of T, and is T_prime = T(p) - Tg(tau,p).




for i = 1:length(tau1000sc);
       similar = find(taurange<=tau1000sc(i)+0.001 & taurange >= tau1000sc(i)-0.001);
       tam = size(similar);
       if(tam(2)>0)
               closest = (tau1000sc(i)-taurange(similar));
               closest = abs(closest);
               closest = (closest == min(closest));
               indclose = find(closest == 1);


               Tgseries(:,i) = SMF_temperature(:,similar(indclose));
       else
               Tgseries(:,i) = NaN.*ones(1,133);
       end;
end;

T=NaN*ones();
for i=1:length(indxf)
  T([indxf(i):indxl(i)])=temp(i);
end


save Tg_series.mat Tgseries
save SMF_temperature.mat SMF_temperature
%save CTD_and_Argo_p_levels.mat dep


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

% PLOTS

% LOAD DATA FROM AOML
read_pies

% THE DATA OLGA PROVIDED ME IS OUTDATED, SO I MUST KNOW IT'S SIZE TO LIMIT PLOTS
%temsize = size(tem);

%--------------
% FOR 2009 ONLY
temsize = [289 NaN];
%--------------



% 1 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
subtightplot(3,1,1,[0.12  0.1],[0.08 0.08], [0.08 0.08])
contourf(timeaxis(1:temsize(1)),-presrange(1:35),Tseries(1:35,1:temsize(1)), 20);
xlabel('Data');
datetick; colorbar;
ylabel('p (dB)');
title('MEU')

% 2 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
subtightplot(3,1,2,[0.12  0.1],[0.08 0.08], [0.08 0.08])


contourf(timeaxis(1:temsize(1)),-presrange(1:35),tem(1:temsize(1),1:35)', 20);

xlabel('Data');
datetick; colorbar;
ylabel('p (dB)');
title('CHRIS')

% 3 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
subtightplot(3,1,3,[0.12  0.1],[0.08 0.08], [0.08 0.08])
contourf(timeaxis(1:temsize(1)),-presrange(1:35),tem(1:temsize(1),1:35)'-Tseries(1:35,1:temsize(1)), 20);


xlabel('Data');
datetick; colorbar;
ylabel('p (dB)');
title('DIF')

print -dpng ds_dif_of_GEMs.png
%!mv ds_dif_of_GEMs.png ~/Dropbox/olga2matheus' (1)'/mestrado/graficos_gem/.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





