path = (['/Volumes/data2/woa13/tem/']);
pathS = (['/Volumes/data2/woa13/sal/']);
more('off');
% NAMES OF PIES MOORINGS
letra = ['ABCD'];

% LOADING DATA FROM EACH MONTH
load pies_latlon;
load Tseries;
load Sseries;
for j = 1:4;
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
end;

T_A = Tseries(:,:,1);
T_B = Tseries(:,:,2);
T_C = Tseries(:,:,3);
T_D = Tseries(:,:,4);
S_D = Sseries(:,:,4);
S_C = Sseries(:,:,3);
S_B = Sseries(:,:,2);
S_A = Sseries(:,:,1);

letra = ['ABCD'];
cor = ['mbgkyr']
figure(1);
for i = 1:4;
	str = (['length_s = size(S_' letra(i) ');']);
	eval(str)
	%subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
	subtightplot(2,2,i,0.1,[0.1 0.05],[0.1 0.05]);
	hold on;
	str2 = (['plot(S_' letra(i) '(:),T_' letra(i) '(:),''.k'',''linewidth'',4);']);
	eval(str2);
	str=(['plot(S_woa_PIES' letra(i) '(:),T_woa_PIES' letra(i) '(:),''.r'',''linewidth'',2);']);
	eval(str)
	if (i == 1 || i == 3)
		ylabel('T (^oC)');
	elseif(i == 3 || i == 4 )
		xlabel('S')
	end;
	str = (['title(''PIES - ' letra(i) ''')']);
	eval(str)
end;

print -depsc2 -r150 woaTS.eps

