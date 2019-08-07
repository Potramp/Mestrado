%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATES POTENTIAL VORTICITY FROM PV_in_vars
% DETECTS MODE WATER USING 1.5xor2e-10 AS LIMIT.
%
%
%
%
% CORTEZI 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; 
clear all;
more('off');

%%% LOADING VARIABLES S, T AND Pdens
load PV_in_vars;
%load PV_in_vars;
load pies_latlon;

%%% CALCULATING VERTICAL DENSITY GRADIENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Vertical gradient of rho_theta... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load CTD_and_Argo_p_levels;

%%%% CUT VARS BELOW 800 dbar
Pdenseries = Pdenseries(1:41,:,:);
Tseries = Tseries(1:41,:,:);
Sseries = Sseries(1:41,:,:);
dep = dep(1:41);

%%%% CUT VARS BELOW 800 dbar
Pdenseries = Pdenseries(1:41,:,:);
Tseries = Tseries(1:41,:,:);
Sseries = Sseries(1:41,:,:);
dep = dep(1:41);

%%%% SMOOTHING T AND Pdens
for j = 1:4
        for i = 1:length(dep);
                sTseries(i,:,j) = conv(Tseries(i,:,j),blackman(31)/sum(blackman(31)),'same');
                sPdenseries(i,:,j) = conv(Pdenseries(i,:,j),blackman(31)/sum(blackman(31)),'same');
        end;
end;

Tseries = sTseries;
Pdenseries = sPdenseries;

clear sTseries sPdenseries;

%%%% making a matrix the same size as the temporal series but for depths(pressure)
zz = repmat(dep',[1 2429 4]);

%%%% making depth out of pressure
zz = sw_dpth(zz,pies_lat_lon(1,1));

%%%% zz is depth, so should be negative.
rho_vgrad = diff(Pdenseries,1,1)./diff(-zz);
%f = f0 + beta.*y;

%%% CALCULATING F (CORIOLIS PARAMETER)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Coriolis parameter... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = coriolisf(pies_lat_lon(1,1));

%%% CALCULATING PV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Potential Vorticity... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% rho_vgrad IS CALCULATED BETWEEN dep LAYERS.
%%% A MEAN DENSITY MUST BE CALCULATED IN ORDER
%%% TO OBTAIN POTENTIAL VORTICITY.

mean_rho = (Pdenseries(1:end-1,:,:) + Pdenseries(2:end,:,:))./2;
PV = f./mean_rho.*rho_vgrad;

%%% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('5. PLots... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
letra = ['ABCD'];
%%%% PV

%%%% MAKE ALL 4 STATIONS WITH FOR
%for j = 1:4;
%	figure(j)
%	subplot(3,1,1);
%	mean_z = (zz(1:end-1,:,:) + zz(2:end,:,:))./2;
%	[cv ch] = contourf(timeaxis,-mean_z(:,1,1),PV(:,:,j),25); colorbar; datetick; ylim([-500 0]);
%	set(ch,'edgecolor','none');
%	%imagesc(timeaxis,mean_z(:,1,1),PV(:,:,j)); colorbar; datetick; ylim([0 500 ]);
%	title('pv')
%
%	%%%% POTENTIAL DENSITY
%	subplot(3,1,2);
%	[cv ch] = contourf(timeaxis,-dep,Pdenseries(:,:,j),25); colorbar; datetick; ylim([-500 0]);
%	set(ch,'edgecolor','none');
%	title('pot. dens.')
%	caxis([1024.6 1027.3])
%	%%%% TEMPERATURE
%	subplot(3,1,3);
%	[cv ch] = contourf(timeaxis,-dep,Tseries(:,:,j),25); colorbar; datetick; ylim([-500 0]);
%	set(ch,'edgecolor','none');
%	hold on;
%	contour(timeaxis,-dep,Tseries(:,:,j),[11.5 11.5],'linecolor','k');        
%	contour(timeaxis,-dep,Tseries(:,:,j),[18.5 18.5],'linecolor','k');        
%	contour(timeaxis,-dep,Tseries(:,:,j),[13 13],'linecolor','m');        
%	contour(timeaxis,-dep,Tseries(:,:,j),[16 16],'linecolor','m');        
%	title('temp')
%	str = (['print -dpng PV_Pden_T_' letra(j) '.png;']);
%	eval(str);
%end;


%%% Detect mode water. Pv limit = pvcrit;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['6. Detecting mode water with pv limit = pvcrit... '])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% STARTING VARIABLES TO SAVE COMPUTING TIME
%stmw = NaN.*ones(size(PV));
%stmw_temp = NaN.*ones(size(PV));
stmw_PV = NaN.*ones(size(PV));

%%%% PV CRITERIUM 
%stmw_pv = (PV <= 1.5e-10);
pvcrit = 1.53e-10;
stmw_pv = (PV <= pvcrit);

%%%% TO PLOT TEMP OF STMW THAT ATTENDS TO BOTH CRITERIA
mean_Tseries = (Tseries(2:end,:,:) + Tseries(1:end-1,:,:))./2;

%%%% TEMP CRITERIUM 
%stmw_temp = (mean_Tseries >= 11.5 & mean_Tseries <=18.5);
%stmw_temp = (mean_Tseries >= 10);


%%%% PV CRITERIUM + TEMP CRITERIUM
%stmw = (stmw_pv == 1 & stmw_pv == stmw_temp);



% Trying to define mode water in terms of pv

%%%%%%%%%%%%%%%%%%%% LOOK AT THIS WHEN CHOOSING PV CRITERIUM %%%%%%%%%%%%%%%%%%
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

%%%% THIS IS FOR ALL CRITERIA
%stmw_PV = PV.*stmw;

%%%% THIS IS WHEN YOU WANT ONLY PV CRITERIA
stmw_PV = PV.*stmw_pv;
stmw_PV(stmw_PV == 0) = NaN;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AND IN TERMS OF TEMP
%stmw_temp = mean_Tseries.*stmw;

% OR USING ONLY PV
stmw_temp = mean_Tseries.*stmw_pv;

	%%%% THIS WILL TAKE ZEROES OUT OF PLOT
	zer = find(stmw_temp == 0);
	stmw_temp(zer) = NaN;


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. PLOTS -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

%for j = 1:4;
%	figure(5);
%	subplot(4,1,j);
%	hold on;
%	contourf(timeaxis,-mean_z(:,1,1),stmw_PV(:,:,j));%,[0e-10 2e-10], 'linecolor', 'k')
%	colorbar; datetick; ylim([-500 0]);
%	caxis([0e-10 2e-10 ]);
%	%contour(timeaxis,-mean_z(:,1,1),stmw_PV(:,:,j),[1.5e-10 1.5e-10], 'linecolor', 'k');colorbar; datetick; ylim([-500 0]);
%	%[cv ch] = contourf(timeaxis,-mean_z(:,1,1),PV(:,:,j),25); colorbar; datetick; ylim([-500 0]);
%
%	figure(6);
%	subplot(4,1,j);
%	% Trying to define mode water in terms of T
%	contourf(timeaxis,-mean_z(:,1,1),stmw_temp(:,:,j));%,[1 1], 'linecolor', 'k');colorbar; datetick; ylim([-500 0]);
%	colorbar; datetick; ylim([-500 0]);
%end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('7. Limited T, and PV values...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% T
stmw_18 = (mean_Tseries >=13 & mean_Tseries <= 19);
t18 = mean_Tseries.*stmw_18;
t18(t18 == 0) = NaN;


%%%% PV
pv18 = PV.*stmw_18;
pv18(pv18 == 0) = NaN;



%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. PLOTS -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%% Temperatures limited by stmw_18
for i =1:4
	figure(7)
	subplot(4,1,i);
	contourf(timeaxis,-mean_z(:,1,1),t18(:,:,i),6); 
	colorbar; datetick('x','yyyy'); ylim([-500 0]);
	title(['T <= 19 and T >= 13 :- PIES ' letra(i)])
end; 

%%%% PV of the points that agree with stmw_18
for i =1:4
	figure(8)
	subplot(4,1,i);
	contourf(timeaxis,-mean_z(:,1,1),pv18(:,:,i),30); 
	colorbar; datetick('x','yyyy'); ylim([-500 0]);
	title(['PV of water with T criterium :- PIES ' letra(i)])
	caxis([0 2e-10])

end; 

%%%% PV and T criteria on T
tpv18 = mean_Tseries.*stmw_18.*stmw_pv;
tpv18(tpv18 == 0) = NaN;
for i =1:4
        figure(9)
        subplot(4,1,i);
        contourf(timeaxis,-mean_z(:,1,1),tpv18(:,:,i),10);
        colorbar; datetick('x','yyyy'); ylim([-500 0]);
        title(['PIES ' letra(i)])
        %%%%%%%%%caxis([0 1.4e-10])

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('8. Calculating vertical gradient of T...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CALCULATING POTENTIAL TEMP.
for i = 1:4
        potTseries(:,:,i) = sw_ptmp(Sseries(:,:,i),Tseries(:,:,i),dep',0);
end;

%%%% making a matrix the same size as the temporal series but for depths(pressure)
zz = repmat(dep',[1 2429 4]);

%%%% making depth out of pressure
zz = sw_dpth(zz,pies_lat_lon(1,1));

%%%% zz is depth, so should be negative.
T_vgrad = diff(potTseries,1,1)./diff(-zz);

%%%% T_VGRAD CRITERIUM
stmw_tgrad = (T_vgrad <= 0.025 );
stmw_tgradpv = PV.*stmw_tgrad;
stmw_tgradpv(stmw_tgradpv == 0) = NaN;

%%%% PV OF WATER THAT FITS TGRAD CRIT
for i =1:4
        figure(10)
        subplot(4,1,i);
        contourf(timeaxis,-mean_z(:,1,1),stmw_tgradpv(:,:,i),35);
        colorbar; datetick('x','yyyy'); ylim([-500 0]);
        title(['Tgrad and pv crit :- PIES ' letra(i)])
        caxis([0 2e-10])

end;

%%%% T OF WATER THAT FITS TGRAD CRIT
mpT = (potTseries(1:end-1,:,:) + potTseries(2:end,:,:))./2; 
mpT_tg = mpT.*stmw_tgrad;
mpT_tg(mpT_tg == 0) = NaN;
for i =1:4
        figure(11)
        subplot(4,1,i);
        contourf(timeaxis,-mean_z(:,1,1),mpT_tg(:,:,i),35);
        colorbar; datetick('x','yyyy'); ylim([-500 0]);
        title(['Tgrad crit on T :- PIES ' letra(i)])
        %caxis([0 2e-10])

end;

%%%% water that fits tgrad and pv crit;
pvc = (PV <= 1.5e-10);
mpT_pv_tg = mpT.*pvc.*stmw_tgrad;
mpT_pv_tg(mpT_pv_tg == 0) = NaN;

for i =1:4
        figure(12)
        subplot(4,1,i);
        contourf(timeaxis,-mean_z(:,1,1),mpT_pv_tg(:,:,i),35);
        colorbar; datetick('x','yyyy'); ylim([-500 0]);
        title(['Tgrad and PV crit on T :- PIES ' letra(i)])
        %caxis([0 2e-10])

end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
