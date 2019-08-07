%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND THE CLOSEST POINT ON ISAS DATA TO EACH PIES 
% STATION. WE THINK THE 16C ISOTHERM LOOKS ODD 
% ON PIES DATA. THIS WILL GIVE AN AVERAGE DEPTH
% FOR THE ISOTHERMAL AS WELL AS A TIME SERIES.
%
% ALSO A TIME SERIES OF T ON THE PIES SITES.
%
%
% CORTEZI, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
more('off');

% DIRECTORY CONTAINING ISAS DATA
path = (['~/Desktop/piero_dados/']);

% NAMES OF PIES MOORINGS
letra = ['ABCD'];

% LOADING DATA FROM EACH MONTH
load pies_latlon;
%for j = 1:4;

% DEFINING DEPTH FROM WHAT PIERO TOLD ME
dep = 0:5:1000;

doitall = input('Run all code?  ');
if( doitall == 1);

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Coordinates for pies points on ISAS data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOAD ANY FILE, JUST FOR COORDINATES;
load([path,'200201']);

% FIND THE NEAREST POINT AVAILABLE IN ISAS TO THE (j) PIES MOORING
for j = 1:4;
                latp = pies_lat_lon(j,1);
                lonp = pies_lat_lon(j,2);
                di = abs(latp-lat);
                closest = find(di == min(di));
                %latp(j) = closest(1);
		indlat(j) = closest;                
		di = abs(lonp-lon);
                closest = find(di == min(di));
		indlon(j) = closest;                
                %lonp(j) = closest(1);
end;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Time Series of T...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% ISAS DATA IS PTEM(LON,LAT,DEPTH)
year = [2002:2014];

%%%% TO TURN SUM TO AVERAGE
cont = 0.*ones(4,12);

%%%% START VARIABLE TO SAVE TIME
month_sum = 0.*ones(201,12,4);

%%%% TIME SERIES FOR ISAS DATA ON MY PIES COORDINATES
ISAS_Tseries = NaN.*ones(201,length(year)*12,4);
depth6 = NaN.*ones(156,4);

%%%% EACH PIES
for k = 2:4;

	%%%% DISPLAY CURRENT PIES
	disp(['PIES ', letra(k), ' !!!']);

	contad = 0;
	%%%% EACH YEAR
	for j = 1:length(year);
		
		%%%% EACH MONTH
		for i = 1:12;
			if(i<10)
				load([path,num2str(year(j)),'0',num2str(i)]);
			else;	
				load([path,num2str(year(j)),num2str(i)]);
			end;

			%%%% TIME SERIES
			junk = squeeze(ptem(indlon(k),indlat(k),:));
			depth16(i+contad,k) = interp1(junk,-dep,16);
			ISAS_Tseries(:,i+contad,k) = junk;
			clear junk

			%%%% MONTH AVG
			month_sum(:,i,k) = month_sum(:,i,k) +squeeze(ptem(indlon(k),indlat(k),:));
			cont(k,i) = cont(k,i)+1;
		end;
		contad = contad + 12;
	end;

	%%%% turn sum into average
	for i = 1:12
		month_avg(:,i,k) = month_sum(:,i,k)./cont(k,i);
	end;
		
end;


%-.-.-.-.-.-.-.-.-. PAUSE -.-.-.-.-.-.-.-.
disp('A Pause before the plots... Press any jey to proceed.')
%pause;

end;


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Plot...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINING DEPTH FROM WHAT PIERO TOLD ME
dep = 0:5:1000;


%%%% FIRST PIES FALLS ON LAND
for i = 1:3;
	subplot(3,1,i)
	contourf([1:12],-dep,month_avg(:,:,i+1),25);colorbar;
	hold on;
	contour([1:12],-dep,month_avg(:,:,i+1),[16 16],'linecolor','m');%colorbar;
end;

%%%% SAVE FIGURE AND MOVE TO OLGA'S FOLDER
print -dpng piero_ISAS.png
%!mv piero_ISAS.png ~/Dropbox/olga2matheus' (1)'/mestrado/graficos_gem/.

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. My data climatology... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load PV_in_vars
month = str2num(datestr(timeaxis,'mm'));
Tsize = size(Tseries);
Tmean = 0.*ones(Tsize(1), 12, Tsize(3));

for j = 1:4;
	for i = 1:12
		indmon = find(month == i);
		Tmean(:,i,j) = nanmean(Tseries(:,indmon,j),2);   
	end;
end;


figure;

%%%% PLOT IT
for i = 1:4;
        subplot(4,1,i)
        contourf([1:12],-dep,Tmean(:,:,i),25);colorbar;
        hold on;
        contour([1:12],-dep,Tmean(:,:,i),[16 16],'linecolor','m');%colorbar;
	ylim([-1000 0]);
end;

%%%% SAVE FIGURE AND MOVE TO OLGA'S FOLDER
print -dpng 16degisotherm.png
%!mv 16degisotherm.png ~/Dropbox/olga2matheus' (1)'/mestrado/graficos_gem/.


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('5. Time series of Pieros data... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETING DATE VECTOR
cont = 1;
for i=1:13;
	for j = 1:12;
		pdate(cont,:) = [year(i) j];
		cont = cont+1;
	end;
end;
pdate = num2str(pdate); pdate = datenum(pdate);  

% DEFINING DEPTH FROM WHAT PIERO TOLD ME
dep = 0:5:1000;

figure;
for i = 2:4;
	subplot(4,1,i)
	%%%% time series
	contourf(pdate,-dep,ISAS_Tseries(:,:,i),25);colorbar;
	hold on;
	%%%% 16deg isotherm
	contour(pdate,-dep,ISAS_Tseries(:,:,i),[16 16],'linecolor','m');
	xlim([pdate(96-12) pdate(96+5*12)]);datetick('x',12,'keeplimits');
end;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('6. Find 16deg isotherm on both data...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%size_piero = size(ISAS_Tseries);
%piero16 = NaN.*ones(size_piero);
%ind16 = find(ISAS_Tseries <= 16.005 & ISAS_Tseries >= 15.995);

%%%% creating same size depth variable to find depth of isotherm;
%zz = repmat(dep',[1 156 4]);

%depth16 = zz(ind16);

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

