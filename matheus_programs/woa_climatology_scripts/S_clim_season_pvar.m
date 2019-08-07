%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This calculates the percent variance of annual and seminannual cicles (of S?)
% for each level of pressupre (p - db) of every PIES from SAM.
%
% 
%
%
%
% CORTEZI, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
more('off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Path to data...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DIRECTORY CONTAINING WOA13 DATA
path = (['/Volumes/data2/woa13/sal/']);
%d = dir([path, '*.mat']);

% NAMES OF PIES MOORINGS
letra = ['ABCD'];

% LOADING DATA FROM EACH MONTH
load pies_latlon;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Fetching points close to PIES ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIND THE NEAREST POINT AVAILABLE IN WOA13 TO THE (j) PIES MOORING
% AND COMPOSING THE AVERAGE YEAR
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
t = 1:12;
table = NaN.*ones(57,2,4);
amplitude_sinus = NaN.*ones(4,57);
for j = 1:4;
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

	
	%%%% For each level of p
	for k = 1:length(dep_woa)
		%%%%%%First T = 6months
		%disp(['PIES-' letra(j) ' depth ' num2str(dep_woa(k)) 'm...']);
		T = 6;
		%%%%% Fit a T-period sine
		str=(['A = sinfitB(1:12,S_woa_PIES' letra(j) '(k,:),T);']);
	        eval(str);
		amplitude_sinus(j,k) = sqrt(A(2)^2 + A(3)^2);
		newx = A(1)+A(2).*t+A(3).*sin((2*pi/T).*t)+A(4).*cos((2*pi/T).*t);
		str = (['pv = pvar(S_woa_PIES' letra(j) '(k,:),newx);']);	
	        eval(str);
		table(k,1,j) = pv;
		%disp([num2str(T) ' months period explains ' num2str(pv) ' percent variance...']);
	
		%%%%%% Second T = 12months
	        T = 12;
	        %%%%% Fit a T-period sine
	        str=(['A = sinfitB(1:12,S_woa_PIES' letra(j) '(k,:),T);']); 
	        eval(str);
	        newx = A(1)+A(2).*t+A(3).*sin((2*pi/T).*t)+A(4).*cos((2*pi/T).*t);      
	        str = (['pv = pvar(S_woa_PIES' letra(j) '(k,:),newx);']);    
	        eval(str);
		table(k,2,j) = pv;
	        %disp([num2str(T) ' months period explains ' num2str(pv) ' percent variance...']);
	end;
end;


cores = 'mkbg';
for i = 1:4
subplot(1,2,1)
title('Semianual - A - m, B - k, C - b, D - g.')
hold on;
str = (['plot(table(:,1,i),-dep_woa, ''' cores(i) ''');']);
eval(str);
hold off;
set(gca,'YTick',[-1500:100:0]);
subplot(1,2,2)
title('Anual')
hold on;
str = (['plot(table(:,2,i),-dep_woa, ''' cores(i) ''');']);
eval(str);
set(gca,'YTick',[-1500:100:0]);
end;



clear i;
figure
cores = 'mkbg';
for i = [1:4];
subplot(1,2,1)
title('AMPLITUDE A - m, B - k, C - b, D - g.')
hold on;
str = (['plot(amplitude_sinus(i,:),-dep_woa, ''' cores(i) ''')']);
disp(['right...the color is ' cores(i)]);
eval(str);
hold off;
set(gca,'YTick',[-1500:100:0]);
subplot(1,2,2)
title('VARIANCE A - m, B - k, C - b, D - g.')
hold on;
str = (['plot(var(S_woa_PIES' letra(i) ',0,2),-dep_woa, ''' cores(i) ''');']);
%str = (['plot(var(S_woa_PIES' letra(i) ',0,2)./mean(S_woa_PIES' letra(i) ',2),-dep_woa, ''' cores(i) ''');']);
eval(str);
set(gca,'YTick',[-1500:100:0]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
