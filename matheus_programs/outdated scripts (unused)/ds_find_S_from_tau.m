% THIS CODE PLOTS A TIME SERIES FOR THE YEAR 2009 FROM THE DATA IN SAM_PIES_data.txt
% USING THE GEM FIELD DONE BY do_temper_field.m

% CREATES GEM
ds_do_sal_field

%LOADS PIES DATA
load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/PIES_data/SAM_PIES_data.txt

% TAU DATA FOR PIES-D
tau_pies = SAM_PIES_data(:,10);

% DATE VECTOR (yy,mm,dd)
date_pies = SAM_PIES_data(:,1:3);
numdate_pies = datenum(date_pies);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- FOR 2009 ONLY TIME SERIES -.-.-.-.-.-.-.-.-.-.-.-
%% SELECTING ONLY THE FIRST YEAR
%ind2009 = find(date_pies(:,1) == 2009);
%tau2009 = tau_pies(ind2009);

%% LOADING VARIABLE TO SAVE COMPUTING TIME
%Sseries = NaN.*ones(133,length(tau2009));

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
%		Sseries(:,i) = SMF_salinity(:,similar(indclose));
%	else
%		Sseries(:,i) = NaN.*ones(1,133);
%	end;
%end;
%timeaxis = numdate_pies(1:289)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- FOR ALL TIME SERIES -.-.-.-.-.-.-.-.-.-.-.-

Ssize = size(tau_pies);

% LOADING VARIABLE TO SAVE COMPUTING TIME
Sseries = NaN.*ones(133,Ssize(1));

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
               Sseries(:,i) = SMF_salinity(:,similar(indclose));
       else
               Sseries(:,i) = NaN.*ones(1,133);
       end;
end;

timeaxis = numdate_pies;



%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
% PLOTS

% LOAD DATA FROM AOML
read_pies

% THE DATA OLGA PROVIDED ME IS OUTDATED, SO I MUST KNOW IT'S SIZE TO LIMIT PLOTS
%salsize = size(tem);

%--------------
% FOR 2009 ONLY
salsize = [289 NaN];
%--------------



% 1 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
subtightplot(3,1,1,[0.12  0.1],[0.08 0.08], [0.08 0.08])
contourf(timeaxis(1:salsize(1)),-presrange(1:35),Sseries(1:35,1:salsize(1)), 20);
xlabel('Data');
datetick; colorbar;
ylabel('p (dB)');
title('MEU')

% 2 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
subtightplot(3,1,2,[0.12  0.1],[0.08 0.08], [0.08 0.08])


contourf(timeaxis(1:salsize(1)),-presrange(1:35),tem(1:salsize(1),1:35)', 20);

xlabel('Data');
datetick; colorbar;
ylabel('p (dB)');
title('CHRIS')

% 3 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
subtightplot(3,1,3,[0.12  0.1],[0.08 0.08], [0.08 0.08])
contourf(timeaxis(1:salsize(1)),-presrange(1:35),tem(1:salsize(1),1:35)'-Sseries(1:35,1:salsize(1)), 20);


xlabel('Data');
datetick; colorbar;
ylabel('p (dB)');
title('DIF')

print -depsc dif_of_GEMs.eps
%!mv dif_of_GEMs.eps ~/Dropbox/olga2matheus' (1)'/mestrado/graficos_gem/.
