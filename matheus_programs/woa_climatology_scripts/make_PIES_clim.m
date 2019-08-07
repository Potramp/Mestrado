%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET CLIMATOLOGIE FROM WOA13 FROM ORION ALREADY MOUNTED
% THIS WILL RESULT A CLIMATOLOGIC YEAR FOR THE  ~ALL~ PIES ECHO SOUNDER
%
% edit: calculates average permanent thermocline and seasonal thermocline profiles.
%	This is used with auto_PIES.m, to determine the temperature gradient we should
%	use to detect mode water.
% - CORTEZI, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;


% DIRECTORY CONTAINING WOA13 DATA
path = (['/Volumes/data2/woa13/tem/']);
pathS = (['/Volumes/data2/woa13/sal/']);
%d = dir([path, '*.mat']);

% NAMES OF PIES MOORINGS
letra = ['ABCD'];

% LOADING DATA FROM EACH MONTH
load pies_latlon;
for j = 1:4;

% FIND THE NEAREST POINT AVAILABLE IN WOA13 TO THE (j) PIES MOORING
% AND COMPOSING THE AVERAGE YEAR
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-		
	for i = 1:12
		if (i<10)
			load([path, 'woa13_tan0' num2str(i) '_04.mat']);
			load([pathS, 'woa13_tan0' num2str(i) '_04.mat']);
		else
			load([path, 'woa13_tan' num2str(i) '_04.mat']);
			load([pathS, 'woa13_tan' num2str(i) '_04.mat']);
		end;

		% FINDING CLOSEST COORDINATES ON WOA FROM THE PIES ECHO SOUNDER j
		latp = pies_lat_lon(j,1);
		lonp = pies_lat_lon(j,2);
		di = abs(latp-lat);
		closest = find(di == min(di));
		latp = closest(1);
		di = abs(lonp-lon);
		closest = find(di == min(di));
		lonp = closest(1);
	
		% LOADING T VALUES FROM T_woa_PIES(j).asc
		% tem = tem(lon,lat,dep)
		str=(['T_woa_PIES' letra(j) '(:,i) = tem(lonp,latp,:);']);
		eval(str);
		str=(['S_woa_PIES' letra(j) '(:,i) = sal(lonp,latp,:);']);
                eval(str);
	end;
	dep_woa = dep;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-		


	%-.-.-.-.-.-.-.-.-.-.-.-.-.-. PLOTS .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
	%
	figure;
	%str=(['surf(1:12,-dep_woa(1:37),T_woa_PIES' letra(j) '(1:37,:),20);']);
	str=(['contourf(1:12,-dep_woa(1:37),T_woa_PIES' letra(j) '(1:37,:),20);']);
	eval(str);
	colorbar;	
	xlabel('month');
	ylabel('depth');
	str = (['title(''PIES ' letra(j) ''');']);
	eval(str);

	% PRINT FIGURE	
	str=(['print -depsc woa13_PIES_' letra(j) '.eps']);
	eval(str);


	clear tem yy mm dd;


% AVERAGE YEAR CALCULATED FROM THE PIES DATA AVAILABLE IN AOML WEBSITE
% EACH MONTH IS SELECTED AND AVERAGED
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-		
	% LOADING T DATA FROM AOML WEBSITE CHANGED BY OLGA
	[yy mm dd jdi dep tem] = read_pie(j);

	
	% DEFINING VARIABLE TO SAVE COMPUTING TIME
	Tmonth= NaN.*ones(133,12);
	for i = 1:12;
		
		% SELECTING EACH MONTH
		month = find(mm == i);
	
		% Tmonth RECEIVES THE MEAN T-PROFILE OF THE i MONTH
		Tmonth(:,i) = nanmean(tem(month,:));
	end;
	%-.-.-.-.-.-.-.-.-.-.-.-.-.-. PLOTS .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
	%
	figure;
	contourf(1:12,-dep(1:35),Tmonth(1:35,:),20);
	colorbar;
	xlabel('month');
	ylabel('depth');
	str = (['title(''PIES ' letra(j) ''');']);
	eval(str);
	% PRINT FIGURE	
	str=(['print -depsc month_mean_PIES_' letra(j) '.eps']);
	eval(str);
end;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Average seasonal and permanent thermocline...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% CALCULATE THE TWO HALF YEAR REGIMES.
%%%% PLOT THE RESULTS
%%%% CALCULTE FIRST DERIVATIVE
for i = 1:4;
	
	%%%% VERTICAL PROFILES OF T
	str = ([letra(i) '_perm = mean(T_woa_PIES' letra(i) '(:,7:12),2);']);
	eval(str);
	str = ([letra(i) '_seaz = mean(T_woa_PIES' letra(i) '(:,1:6),2);']);
	eval(str);
		%% PLOTS
		figure(1)
		subplot(2,2,i)
		str = (['plot(' letra(i) '_seaz, -dep_woa,''linewidth'',3);']);
		eval(str);
		grid on;	
		figure(2)
		subplot(2,2,i);
		str = (['plot(' letra(i) '_perm, -dep_woa,''linewidth'',3);']);
		eval(str);
		grid on;	

	%%% Derivatives
	str = (['deltaT_' letra(i) '_seaz = diff(' letra(i) '_seaz)./diff(-dep_woa)']);
	eval(str);
	str = (['deltaT_' letra(i) '_perm = diff(' letra(i) '_perm)./diff(-dep_woa)']);
	eval(str);
		%% PLOTS
		figure(3)
		subplot(2,2,i)
		str = (['plot(deltaT_' letra(i) '_seaz,-(dep_woa(1:end-1)+dep_woa(2:end))./2);']);
		eval(str);
		grid on;	
		figure(4)
		subplot(2,2,i)
		str = (['plot(deltaT_' letra(i) '_perm,-(dep_woa(1:end-1)+dep_woa(2:end))./2,''linewidth'',3);']);
		eval(str);
		grid on;	
	
end;




