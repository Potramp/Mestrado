% GET CLIMATOLOGI FROM WOA13 FROM ORION ALREADY MOUNTED
% THIS WILL RESULT A CLIMATOLOGIC YEAR FOR THE PIESD ECHO SOUNDER
close all;
clear all;


% DIRECTORY CONTAINING WOA13 DATA
path = (['/Volumes/data2/woa13/tem/']);
%d = dir([path, '*.mat']);


% LOADING DATA FROM EACH MONTH
for i = 1:12
	if (i<10)
		load([path, 'woa13_tan0' num2str(i) '_04.mat']);
	else
		load([path, 'woa13_tan' num2str(i) '_04.mat']);
	end;
	

	% tem = tem(lon,lat,dep)
	T_woa_PIESD(:,i) = tem(102,62,:);
end;
dep_woa = dep;
% LOADING T DATA FROM AOML WEBSITE CHANGED BY OLGA
read_pies

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
figure;
contourf(1:12,-dep_woa(1:37),T_woa_PIESD(1:37,:),20);
