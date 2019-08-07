%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATES POTENTIAL VORTICITY FROM PV_in_vars
% DETECTS MODE WATER USING pvc AS LIMIT.
%
% -automated version for tests of parameters
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

%%%% CUT VARS BELOW 600 dbar
Pdenseries = Pdenseries(1:36,:,:);
Ts = Ts(1:36,:,:);
Ss = Ss(1:36,:,:);
dep = dep(1:36);

%%%% SMOOTHING T AND Pdens
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- BLACKMAN FILTER -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%for j = 1:4
%        for i = 1:length(dep);
%                sTs(i,:,j) = conv(Ts(i,:,j),blackman(31)/sum(blackman(31)),'same');
%                sPdenseries(i,:,j) = conv(Pdenseries(i,:,j),blackman(31)/sum(blackman(31)),'same');
%        end;
%end;

%Ts = sTs;
%Pdenseries = sPdenseries;
%clear sTs sPdenseries;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- /BLACKMAN FILTER -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-


%%%% VERTICAL GRADIENT OF DENSITY
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%% making a matrix the same size as the temporal series but for depths(pressure)
zz = repmat(dep',[1 2429 4]);

%%%% making depth out of pressure
zz = sw_dpth(zz,pies_lat_lon(1,1));

%%%% zz is depth, so should be negative.
rho_vgrad = diff(Pdenseries,1,1)./diff(-zz);
%f = f0 + beta.*y;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. Variables for plots...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% VARIABLES FOR PLOTS
letra = ['ABCD'];
mean_z = (zz(1:end-1,:,:) + zz(2:end,:,:))./2;

%%%% TO PLOT TEMP OF STMW THAT ATTENDS TO BOTH CRITERIA
mean_Ts = (Ts(2:end,:,:) + Ts(1:end-1,:,:))./2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('5. Calculating vertical gradient of T...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CALCULATING POTENTIAL TEMP.
for i = 1:4
        potTs(:,:,i) = sw_ptmp(Ss(:,:,i),Ts(:,:,i),dep',0);
end;

%%%% making a matrix the same size as the temporal series but for depths(pressure)
zz = repmat(dep',[1 2429 4]);

%%%% making depth out of pressure
zz = sw_dpth(zz,pies_lat_lon(1,1));

%%%% zz is depth, so should be negative.
T_vgrad = diff(potTs,1,1)./diff(-zz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('6. STMW criteria... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- USER INPUT -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%% T_VGRAD CRITERIUM
stmw_tgrad = (T_vgrad <= 0.015 );
stmw_tgradpv = PV.*stmw_tgrad;
stmw_tgradpv(stmw_tgradpv == 0) = NaN;

%%%% T OF WATER THAT FITS TGRAD CRIT
mpT = (potTs(1:end-1,:,:) + potTs(2:end,:,:))./2; 
mpT_tg = mpT.*stmw_tgrad;
mpT_tg(mpT_tg == 0) = NaN;

%%%% T LIMITS
minT = [11.5 13]; 
maxT = [19.5 16];

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%% PV CRITERIUM
%%%% water that fits tgrad and pv crit;
%ui = input(['Pv crit (times e-10): ']);
%str = (['pvc(jj).jj)= (PV <= ' num2str(ui) 'e-10);']);
%eval(str);

pvc = [1.5e-10 2e-10];
%pvc = [1:0.1:1.5].*10^(-10);
%%%%% THIS MUST COME INSIDE THE LOOP
%mpT_pv_tg = mpT.*pvc(jj).*stmw_tgrad;
%mpT_pv_tg(mpT_pv_tg == 0) = NaN;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

%%%% ITERARION TO TEST ALL T LIMITS AND PVC
%-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

%%%% LOOPS MINT, MAXT AND PVC
for ii = 1:2
for jj = 1:length(pvc)

%%%% correction needed because of loop
pvboul = (PV <= pvc(jj));
mpT_pv_tg = mpT.*pvboul.*stmw_tgrad;
mpT_pv_tg(mpT_pv_tg == 0) = NaN;

wT = mpT;

%%%% CREATES A DIRECTORY IF IT DOES NOT ALREADY EXISTS
%%%% FIXING PROBLEM WITH STR
junk = num2str(pvc(jj)/1e-11);
junk2 = num2str(floor(minT(ii)));
junk3 = num2str(floor(maxT(ii)));

path = (['~/Dropbox/Mestrado/matheus_programs/without_correction/' junk2 '_' junk3 '_', junk]);
if(exist(path)==0)
	eval(['!mkdir ' path]);
end;

	critT = (mean_Ts >= minT(ii) & mean_Ts <= maxT(ii));
	tcritT = mean_Ts.*critT;
	tcritT(tcritT == 0) = NaN;
	%%%% PV OF WATER UNDER T CRIT
	pvcritT = PV.*critT;
	pvcritT(pvcritT == 0) = NaN;

	%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- /USER INPUT -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


	%%% PLOTS
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp('7. PLots... ')
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%% ON FIGURE 1 I WANT 1 ISOLINE FOR EACH DEGREE
	c = length(minT(ii):1:maxT(ii));

	%%%% Temperatures limited by critT
	for i =1:4
		figure(1)
		subplot(4,1,i);
		contourf(timeaxis,-mean_z(:,1,1),tcritT(:,:,i),c); 
		colorbar; datetick('x','yyyy'); ylim([-600 0]);
		str = (['title([''T <= '' num2str(maxT(ii)) '' and T >= '' num2str(minT(ii)) '' :- PIES '' letra(i)])']);
		eval(str);
	end; 
	str = (['print -dpsc ' path '/fig1.eps']);
	eval(str);

	%%%% PV of the points that agree with critT
	for i =1:4
		figure(2)
		subplot(4,1,i);
		contourf(timeaxis,-mean_z(:,1,1),pvcritT(:,:,i),30); 
		colorbar; datetick('x','yyyy'); ylim([-600 0]);
		title(['PV of water with T criterium :- PIES ' letra(i)])
		caxis([0 2e-10])

	end; 
	str = (['print -dpsc ' path '/fig2.eps']);
	eval(str);

	%%%% PV and T criteria on T
	tpvcritT = mean_Ts.*critT.*pvboul;
	tpvcritT(tpvcritT == 0) = NaN;
	for i =1:4
	        figure(3)
	        subplot(4,1,i);
	        contourf(timeaxis,-mean_z(:,1,1),tpvcritT(:,:,i),10);
	        colorbar; datetick('x','yyyy'); ylim([-600 0]);
	        title(['PIES ' letra(i)])
	        %%%%%%%%%caxis([0 1.4e-10])

	end;
	str = (['print -dpsc ' path '/fig3.eps']);
	eval(str);

	%%%% PV OF WATER THAT FITS TGRAD CRIT
	for i =1:4
	        figure(4)
	        subplot(4,1,i);
	        contourf(timeaxis,-mean_z(:,1,1),stmw_tgradpv(:,:,i),35);
	        colorbar; datetick('x','yyyy'); ylim([-600 0]);
	        title(['Tgrad and pv crit :- PIES ' letra(i)])
	        caxis([0 2e-10])

	end;
	str = (['print -dpsc ' path '/fig4.eps']);
	eval(str);


	%%%% T OF WATER THAT FITS TGRAD CRIT
	for i =1:4
	        figure(5)
	        subplot(4,1,i);
	        contourf(timeaxis,-mean_z(:,1,1),mpT_tg(:,:,i),35);
	        colorbar; datetick('x','yyyy'); ylim([-600 0]);
	        title(['Tgrad crit on T :- PIES ' letra(i)])
	        %caxis([0 2e-10])

	end;
	str = (['print -dpsc ' path '/fig5.eps']);
	eval(str);

	for i =1:4
	        figure(6)
	        subplot(4,1,i);
	        contourf(timeaxis,-mean_z(:,1,1),mpT_pv_tg(:,:,i),35);
	        colorbar; datetick('x','yyyy'); ylim([-600 0]);
	        title(['Tgrad and PV crit on T :- PIES ' letra(i)])
	        caxis([0 15])

	end;
	str = (['print -dpsc ' path '/fig6.eps']);
	eval(str);

	%%%% WATER THAT FITS PV, TGRAD AND T CRITERIA
	w = mpT.*pvboul.*stmw_tgrad.*critT;
	w(w==0) = NaN;
	%cn = input('Number of contours on fig 6: ');
	cn = 5;

	for i = 1:4;
		figure(7)
		subplot(4,1,i);
		contourf(timeaxis,-mean_z(:,1,1),w(:,:,i),cn);
		colorbar; datetick('x','yyyy'); ylim([-600 0]);
		title(['Tgrad, PV and T criteria :- PIES ' letra(i)])
	end;
	str = (['print -dpsc ' path '/fig7.eps']);
	eval(str);
	
	%%%% SAVE w
	str = (['save ' path '/w.mat w mean_z PV wT']);
	eval(str);
	
end;
end;
%-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
