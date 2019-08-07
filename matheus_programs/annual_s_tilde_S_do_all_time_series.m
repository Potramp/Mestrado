% THIS CODE PLOTS A TIME SERIES FOR THE YEAR 2009 FROM THE DATA IN SAM_PIES_data.txt
% USING THE GEM FIELD DONE BY ds_do_sal_field.m

% CREATES GEM
ds_do_sal_field

%LOADS PIES DATA
load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/PIES_data/SAM_PIES_data.txt
%load /Users/apple/Dropbox/Mestrado/chris_programs/PIES_data/SAM_PIES_data.txt

% TAU DATA FOR PIES-D
tau_pies = SAM_PIES_data(:,[4 6 8 10]);

% DATE VECTOR (yy,mm,dd)
date_pies = SAM_PIES_data(:,1:3);
numdate_pies = datenum(date_pies);
mm = date_pies(:,2);
yy = date_pies(:,1);
dd = date_pies(:,3);
% S_tilde_y FOR SEASONAL CORRECTION
load S_tilde_y

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
%%% Making S_tilde_y the size as SMF_salinity filling with NaNs.
%junk = size(SMF_salinity);
%junk = 0.*ones(junk(1),12);
%junk2 = size(S_tilde_y);
%junk(1:junk2,:) = S_tilde_y;
%S_tilde_y = junk;

junk = size(SMF_salinity);
junk2 = size(S_tilde_y);
junk = 0.*ones(junk(1),junk2(2));
junk2 = size(S_tilde_y);
junk(1:junk2,:) = S_tilde_y;
S_tilde_y = junk;

clear junk*; 

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. RECONSTRUCTING TIME SERIES -.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reconstructing time series of S for all PIES moorings ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% THIS LOOP WILL MAKE TIME SERIES FOR ALL 4 PIES
%%%% 
for j = 1:4;

	% FINDS THE SMF_salinity PROFILE THAT HAS THE TAU VALUE CLOSEST TO THAT OF THE
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
			%%% I need to make the seasonal correction for S as in Watts 2001.
	       		indyear = datenum(yy(i) - yy(i), mm(i), dd(i));
			if(indyear <= 365)
                                Sseries(:,i,j) = SMF_salinity(:,similar(indclose)) + S_tilde_y(:,indyear);
                        else;
                                Sseries(:,i,j) = SMF_salinity(:,similar(indclose)) + S_tilde_y(:,365);
                        end;


	       else
	               Sseries(:,i,j) = NaN.*ones(1,133);
	       end;
	end;
end;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. RECONSTRUCTING TIME SERIES -.-.-.-.-.-.-.-.-.-.-.-.-.

timeaxis = numdate_pies;

%%%%-.-.-.-.-.-.-.-.-.-.-.-.-.  Calculating Sg -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%% For S seasonal cycle as in Watts2001, we apply the same method did in ds_do_sal_field
%%%% for tau, but this time for S. Sg is calculated from the tau data derived from CTD.
%%%% S_prime  will be the annual march of S, and is S_prime = S(p) - Sg(tau,p).




for i = 1:length(tau1000sc);
       similar = find(taurange<=tau1000sc(i)+0.001 & taurange >= tau1000sc(i)-0.001);
       tam = size(similar);
       if(tam(2)>0)
               closest = (tau1000sc(i)-taurange(similar));
               closest = abs(closest);
               closest = (closest == min(closest));
               indclose = find(closest == 1);


               Sgseries(:,i) = SMF_salinity(:,similar(indclose));
       else
               Sgseries(:,i) = NaN.*ones(1,133);
       end;
end;

S=NaN*ones(length(indxf));
for i=1:length(indxf)
  S([indxf(i):indxl(i)])=sal(i);
end


save Sg_series.mat Sgseries
save SMF_salinity.mat SMF_salinity
%save CTD_and_Argo_p_levels.mat dep


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

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
%salsize = size(sal);

%--------------
% FOR 2009 ONLY
salsize = [289 NaN];
%--------------
letra = ['ABCD'];
for j = 1:4;
	figure(1)

	% 1 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	subtightplot(4,1,j,[0.12  0.1],[0.08 0.08], [0.08 0.08])
	[cs, h] = contourf(timeaxis(1:salsize(1)),-presrange(1:35),Sseries(1:35,1:salsize(1),j), 20);
	xlabel('Data');
	datetick('x', 28); colorbar;
	ylabel('p (dbar)');
	clabel(cs, h);
	title(['PIES ', num2str(letra(j))]);
	
	figure(2)

	% 2 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	subtightplot(4,1,j,[0.12  0.1],[0.08 0.08], [0.08 0.08])
	plot(timeaxis(1:salsize(1)),sal_p(1:salsize(1),1,j));
%	plot(timeaxis(1:salsize(1)),Sseries(1,1:salsize(1),j));
	xlabel('Data');
	datetick('x', 28); %colorbar;
%	ylabel('T (degrees)');
	ylabel('S (psu)')
	title(['PIES ', num2str(letra(j))]);

	figure(3)

	% 2 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	subtightplot(4,1,j,[0.12  0.1],[0.08 0.08], [0.08 0.08])
%	plot(timeaxis(1:salsize(1)),sal_p(1:salsize(1),1,j));
	plot(timeaxis(1:salsize(1)),Sseries(1,1:salsize(1),j));
	xlabel('Data');
	datetick('x', 28); %colorbar;
%	ylabel('T (degrees)');
	ylabel('S (psu)')
	title(['PIES ', num2str(letra(j))]);


figure(4);
        % 1 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	subtightplot(3,1,1,[0.12  0.1],[0.08 0.08], [0.08 0.08])
        [cs, h] = contourf(timeaxis(1:salsize(1)),-presrange(1:35),Sseries(1:35,1:salsize(1),j), 20);
        xlabel('Data');
        datetick('x', 28); colorbar;
        ylabel('p (dbar)');
%        clabel(cs, h);
        title(['PIES ', num2str(letra(j))]);




	% 2 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	subtightplot(3,1,2,[0.12  0.1],[0.08 0.08], [0.08 0.08])
	
	
	contourf(timeaxis(1:salsize(1)),-presrange(1:35),sal(1:salsize(1),1:35)', 20);
	
	xlabel('Data');
	datetick('x', 28); colorbar;
	ylabel('p (dB)');
	title('CHRIS')
	
	% 3 -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
	subtightplot(3,1,3,[0.12  0.1],[0.08 0.08], [0.08 0.08])
	contourf(timeaxis(1:salsize(1)),-presrange(1:35),sal(1:salsize(1),1:35)'-Sseries(1:35,1:salsize(1)), 20);
	
	
	xlabel('Data');
	datetick('x', 28); colorbar;
	ylabel('p (dbar)');
	title('DIF')
	

end;

save Sseries.mat Sseries;
%print -dpng time_series_PIES.png
%!mv time_series_PIES.png ~/Dropbox/olga2matheus' (1)'/mestrado/graficos_gem/.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





