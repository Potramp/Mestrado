% THIS CODE PLOTS A TIME SERIES FOR THE YEAR 2009 FROM THE DATA IN SAM_PIES_data.txt
% USING THE GEM FIELD DONE BY ds_do_sal_field.m

% CREATES GEM
do_sal_mat

%LOADS PIES DATA
load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/PIES_data/SAM_PIES_data.txt
%load /Users/apple/Dropbox/Mestrado/chris_programs/PIES_data/SAM_PIES_data.txt

% TAU DATA FOR PIES-D
tau_pies = SAM_PIES_data(:,[4 6 8 10]);

% DATE VECTOR (yy,mm,dd)
date_pies = SAM_PIES_data(:,1:3);
numdate_pies = datenum(date_pies);
mm = date_pies(:,2);
% S_tilde FOR SEASONAL CORRECTION
%load S_tilde

%%%% THIS SECTION IS NOT BEING USED ANYMORE. 
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
Sseries = 0.*ones(133,Ssize(1));


%%% PREPARING FOR SEASONAL CORRECTION
%%% Making S_tilde the size as SMF_salinity filling with NaNs.
%junk = size(SMF_salinity);
%junk = 0.*ones(junk(1),12);
%junk2 = size(S_tilde);
%junk(1:junk2,:) = S_tilde;
%S_tilde = junk;

%clear junk*; 

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. RECONSTRUCTING TIME SERIES -.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reconstructing time series of T for all PIES moorings ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% THIS LOOP WILL MAKE TIME SERIES FOR ALL 4 PIES
%%%% 
for j = 1:4;

	% FINDS THE SMF_Temperature PROFILE THAT HAS THE TAU VALUE CLOSEST TO THAT OF THE
	% PIES DATA
	for i = 1:length(tau_pies(:,j));
	       similar = find(taurange<=tau_pies(i,j)+0.001 & taurange >= tau_pies(i,j)-0.001);
	       tam = size(similar);
	       if(tam(2)>0)
	               closest = (tau_pies(i,j)-taurange(similar));
	               closest = abs(closest);
	               closest = (closest == min(closest));
	               indclose = find(closest == 1);
			%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. SEASONAL CORRECTION -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
			%%% I need to make the seasonal correction for T as in Watts 2001.
	               %Sseries(:,i,j) = SMF_salinity(:,similar(indclose)) + S_tilde(:,mm(i));
	               Sseries(:,i,j) = SMF_salinity(:,similar(indclose));
	       else
	               Sseries(:,i,j) = NaN.*ones(1,133);
	       end;
	end;
end;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. RECONSTRUCTING TIME SERIES -.-.-.-.-.-.-.-.-.-.-.-.-.

timeaxis = numdate_pies;
% PLOTS

% LOAD DATA FROM AOML
%%%% READ PIES-D
read_pies
letra = ['ABCD'];
%%%% READ ALL PIES
for j = 1:4;
	str = (['[yy_p mm_p dd_p jd_p dep_p sal_p_' letra(j) '] = read_pie(j);']);
	eval(str);
	str = (['sal_p(:,:,j) = sal_p_' letra(j) ';']);
	eval(str);
end;


% THE DATA OLGA PROVIDED ME IS OUTDATED, SO I MUST KNOW IT'S SIZE TO LIMIT PLOTS
salsize = size(sal);

%--------------
% FOR 2009 ONLY
%salsize = [289 NaN];
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

letra = ['ABCD'];
for j = 1:4;
	figure(1)

	% 1 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	subtightplot(4,1,j,[0.12  0.1],[0.08 0.08], [0.08 0.08])
	contourf(timeaxis(1:salsize(1)),-presrange(1:35),Sseries(1:35,1:salsize(1),j), 20);
	xlabel('Data');
	datetick; colorbar;
	ylabel('p (dbar)');
	title(['PIES ', num2str(letra(j))]);
	
	figure(2)

	% 2 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	subtightplot(4,1,j,[0.12  0.1],[0.08 0.08], [0.08 0.08])
	plot(timeaxis(1:salsize(1)),sal_p(1:salsize(1),1,j));
%	plot(timeaxis(1:salsize(1)),Sseries(1,1:salsize(1),j));
	xlabel('Data');
	datetick; %colorbar;
	ylabel('T (deg)');
	title(['PIES ', num2str(letra(j))]);


	figure(3)

	% 3 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	subtightplot(4,1,j,[0.12  0.1],[0.08 0.08], [0.08 0.08])
%	plot(timeaxis(1:salsize(1)),sal_p(1:salsize(1),1,j));
	plot(timeaxis(1:salsize(1)),Sseries(1,1:salsize(1),j));
	xlabel('Data');
	datetick; %colorbar;
	ylabel('T (deg)');
	title(['PIES ', num2str(letra(j))]);

	save Sseries.mat Sseries;

	

end;
Ss = Sseries(:,:,:);
save Ss.mat Ss
print -dpng lalala.png
%!mv time_series_PIES.png ~/Dropbox/olga2matheus' (1)'/mestrado/graficos_gem/.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





