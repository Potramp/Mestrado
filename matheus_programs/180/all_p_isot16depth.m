%%%% FOR ALL REFERENCE LEVELS

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
%doitall = input('Run all code?  ');
doitall = 1;
if( doitall == 1);

close all;
clear all;
more('off');

% DIRECTORY CONTAINING ISAS DATA
path = (['~/Desktop/piero_dados/']);

% NAMES OF PIES MOORINGS
letra = ['ABCD'];

% LOADING DATA FROM EACH MONTH
load ~/Dropbox/Mestrado/matheus_programs/pies_latlon;
%for j = 1:4;

% DEFINING DEPTH FROM WHAT PIERO TOLD ME
dep_p = 0:5:1000;


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
%depth16 = NaN.*ones(156,4);

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
			%depth16(i+contad,k) = interp1(junk,-dep_p,16);
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

save ~/Dropbox/Mestrado/matheus_programs/ISAS_Tseries.mat ISAS_Tseries;
%-.-.-.-.-.-.-.-.-. PAUSE -.-.-.-.-.-.-.-.
%disp('A Pause before the plots... Press any jey to proceed.')
%pause;

end;


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Plot...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINING DEPTH FROM WHAT PIERO TOLD ME
dep_p = 0:5:1000;


%%%% FIRST PIES FALLS ON LAND
for i = 2:4;
	subplot(4,1,i)
	contourf([1:12],-dep_p,month_avg(:,:,i),25);colorbar;
	hold on;
	contour([1:12],-dep_p,month_avg(:,:,i),[16 16],'linecolor','m');%colorbar;
	title('ISAS data average year');
end;

%%%% SAVE FIGURE AND MOVE TO OLGA'S FOLDER
print -dpng piero_ISAS.png
%!mv piero_ISAS.png ~/Dropbox/olga2matheus' (1)'/mestrado/graficos_gem/.

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. My data climatology... ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ~/Dropbox/Mestrado/matheus_programs/PV_in_vars
load Tseries;
month = str2num(datestr(timeaxis,'mm'));
Tsize = size(Tseries);
Tmean = 0.*ones(Tsize(1), 12, Tsize(3));

rmpath ~/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts
for j = 1:4;
	for i = 1:12
		indmon = find(month == i);
		Tmean(:,i,j) = nanmean(Tseries(:,indmon,j),2);   
	end;
end;



%%%% PLOT IT
figure;

for i = 1:4;
        subplot(4,1,i)
        contourf([1:12],-dep,Tmean(:,:,i),25);colorbar;
        hold on;
        contour([1:12],-dep,Tmean(:,:,i),[16 16],'linecolor','m');%colorbar;
	ylim([-1000 0]);
	title('PIES average year');
end;
hold off;

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
%dep_p = 0:5:1000;

figure;
for i = 2:4;
	subplot(4,1,i)
	%%%% time series
	contourf(pdate,-dep_p,ISAS_Tseries(:,:,i),25);colorbar;
	hold on;
	%%%% 16deg isotherm
	contour(pdate,-dep_p,ISAS_Tseries(:,:,i),[16 16],'linecolor','m');
	xlim([pdate(96-12) pdate(96+5*12)]);datetick('x',12,'keeplimits');
end;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('6. Find 16deg isotherm on both data...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FIXING DEPTH PROBLEM
load ~/Dropbox/Mestrado/matheus_programs/PIES_depth;


%%%% FIND 16 ISOT DEPTH ON MY DATA
zz = repmat(dep',[1 156 4]);

%depth16 = zz(ind16);


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. MY DATA -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
for j = 1:4;
	for i = 1:length(Tseries);
		if(sum(isnan(Tseries(:,i,j)) == 0 ));
			if(sum(diff(Tseries(:,i,j)) ~= 0));
				depth16(j,i) = interp1(Tseries(:,i,j),zz(:,1,1),16);
			else;		
				depth16(j,i) = NaN;
			end;
		else;
			depth16(j,i) = NaN;
		end;
	end;
	
	%%%% -.-.-.-.-.-.-.-.-.-. PLOT -.-.-.-.-.-.-.-.-.-.-.-.
	%figure(1)
	subplot(4,1,j);
	
	[cv ch] = contourf(timeaxis,-dep,Tseries(:,:,j),25);
	set(ch,'edgecolor','none');
	ylim([-500 0]);
	hold on;
	
	%%%% fixing length with Tseries;
	%junk = sum(Tseries(:,:,j),1);
	%indofnan = find(isnan(junk));
	%depth16(indofnan) = NaN; 
		
	
	%figure(2)
	%subplot(4,1,j);
	%plot(timeaxis,depth16(j,:));
end;

hold off;	

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. ISAS DATA -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
% DEFINING DEPTH FROM WHAT PIERO TOLD ME
dep_p = 0:5:1000;
ITsize = size(ISAS_Tseries);
for j = 2:4;
	for i = 1:ITsize(2);
		if(sum(isnan(ISAS_Tseries(:,i,j)) == 0));
			ISAS_depth16(j,i) = interp1(ISAS_Tseries(:,i,j),dep_p,16);
		else;
			ISAS_depth16(j,1) = NaN;
		end;
	end;
	
	%%%% -.-.-.-.-.-.-.-.-.-. PLOT -.-.-.-.-.-.-.-.-.-.-.-.
	subplot(4,1,j);
	
	[cv ch] = contourf(pdate(84:end),-dep_p,ISAS_Tseries(:,84:end,j),25);
	set(ch,'edgecolor','none');
	ylim([-500 0]);
	hold on;

	plot(pdate(84:end),ISAS_depth16(j,84:end));
end;
hold off;

%%%%
save ../ISAS_depth16.mat ISAS_depth16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('7. Dif of 16 deg isothermal, mine and piero''s...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% SET A VARIABLE WITH THE YEARS
m_year = str2num(datestr(timeaxis,'yyyy'));

%%%% STARTING VARIABLE TO SAVE COMPUTING TIME
%monthly_depth16 = NaN.*ones(,156);

%%%% TO CHECK WHEN THE MONTH CHANGES
last_month = month(1);

%cont_m = 1;

%%%% YEAR AND MONTH IN ONE VARIABLE;
m_date = [m_year';month'];
m_date = m_date';

junk_year = 2009:2014;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%monthly_depth = NaN.*ones(4,156);
%dp16size = size(depth16);
num_of_years = length(junk_year);

%%%% THE PIES
for i = 1:4
	%%%% THIS COUNTER WILL PUT THE AVERAGE ON THE RIGHT PLACE
	cont = 1;

	%%%% THE YEAR
	for j = 1:num_of_years;
	
		
		ind_year = find(m_year == junk_year(j));
		%%%% THE MONTH
		for k = 1:12;
			ind_month = find(month(ind_year) == k);
			valid_data = size(ind_month);		

			if(valid_data(1) > 0)
				%%%% TAKE AWAY NAN PROFILES
				if(sum(depth16(i,ind_month)) ~= NaN);
					monthly_depth(i,cont) = nanmean(depth16(i,ind_year(ind_month)));
				else;
					monthly_depth(i,cont) = NaN;
				end;
			else;
				monthly_depth(i,cont) = NaN;
			end;
			
			%%%% UPDATE COUNTER;
			cont = cont+1;

		end;
	end;
end;





%%%% just in case it goes wrong with the nans
%monthly_depth16(:,1) = NaN.*ones(1,length(monthly_depth16));

%%%% DIF OF HIS DATA AND MINE
ISAS_depth16 = ISAS_depth16';

%%%% 85 brcause it must start on January of 2009;
error = ISAS_depth16(85:end,:)'-monthly_depth;

%%%% PLOT OF SAID DIF
%figure(4)
%plot(error(:,2:4))

%%%% SAVING DIF AND FIGURE
save error.mat monthly_depth ISAS_depth16 error depth16;
%print -dpng error.png;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


