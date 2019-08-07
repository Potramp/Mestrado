% GET CLIMATOLOGI FROM WOA13 FROM ORION ALREADY MOUNTED
% THIS WILL RESULT A CLIMATOLOGIC YEAR FOR THE  ~ALL~ PIES ECHO SOUNDER

%%%% !!!! THIS IS FOR SALINITY !!!! 

close all;
clear all;


% DIRECTORY CONTAINING WOA13 DATA
path = (['/Volumes/data2/woa13/sal/']);
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
		else
			load([path, 'woa13_tan' num2str(i) '_04.mat']);
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
	
		% LOADING S VALUES FROM S_woa_PIES(j).asc
		% sal = sal(lon,lat,dep)
		str=(['S_woa_PIES' letra(j) '(:,i) = sal(lonp,latp,:);']);
		eval(str);
	end;
	dep_woa = dep;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-		


	%-.-.-.-.-.-.-.-.-.-.-.-.-.-. PLOTS .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
	%
	figure;
	%str=(['surf(1:12,-dep_woa(1:37),S_woa_PIES' letra(j) '(1:37,:),20);']);
	str=(['contourf(1:12,-dep_woa(1:37),S_woa_PIES' letra(j) '(1:37,:),20);']);
	eval(str);
	colorbar;	
	xlabel('month');
	ylabel('depth');
	%clabel('Salinity [ssu]');
	str = (['title(''PIES ' letra(j) ''');']);
	eval(str);

	% PRINT FIGURE	
	str=(['print -depsc S_woa13_PIES_' letra(j) '.eps']);
	eval(str);


	clear sal yy mm dd;


% AVERAGE YEAR CALCULATED FROM THE PIES DATA AVAILABLE IN AOML WEBSITE
% EACH MONTH IS SELECTED AND AVERAGED
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-		
	% LOADING T DATA FROM AOML WEBSITE CHANGED BY OLGA
	[yy mm dd jd dep sal] = S_read_pie(j);

	
	% DEFINING VARIABLE TO SAVE COMPUTING TIME
	Smonth= NaN.*ones(133,12);
	for i = 1:12;
		
		% SELECTING EACH MONTH
		month = find(mm == i);
	
		% Smonth RECEIVES THE MEAN T-PROFILE OF THE i MONTH
		Smonth(:,i) = nanmean(sal(month,:));
	end;
	%-.-.-.-.-.-.-.-.-.-.-.-.-.-. PLOTS .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
	%
	figure;
	contourf(1:12,-dep(1:35),Smonth(1:35,:),20);
	colorbar;
	xlabel('month');
	ylabel('depth');
	str = (['title(''PIES ' letra(j) ''');']);
	eval(str);
	% PRINT FIGURE	
	str=(['print -depsc S_month_mean_PIES_' letra(j) '.eps']);
	eval(str);
end;
