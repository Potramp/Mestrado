%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w IN THE CURRENT FOLDER CONTAINS ALL THAT IS 
% NEEDED TO DETERMINE MODE WATER THICKNESS. 
% FROM w, PV and mean_z THIS WILL FIND min_d,
% max_d and thick.
%
%
%
%
% 
% -CORTEZI, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% CHECK WHERE YOU ARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Current folder...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pwd

%%%% GET VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Load variables...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load w.mat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. min_d, max_d, thick...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% FIND WHAT WAS JUDGED AS MODE WATER
mask = isnan(w);

%%%% DIMENSIONS OF w
sw = size(w);

%%%% GO THROUGH STATIONS
for i = 1:4;
	%%%% GO THROUGH DAYS
	for j = 1:sw(2);
		%%%% valid data
		v = find(mask(:,j,i) == 0);
		
		%%%% CHECK IF THERE IS MODE WATER AT ALL
		if(length(v) ~= 0 )
			min_d(j,i) = mean_z(v(1),j,i);
			max_d(j,i) = mean_z(v(end),j,i);
		else;
			min_d(j,i) = NaN;
			max_d(j,i) = NaN;
		end;

	end;
end;

thick = max_d-min_d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. PLots')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% STATION NAMES
letra = ['ABCD'];

figure(1)
for i = 1:4
subtightplot(4,1,i,0.08,0.05,0.05);
hold on;
plot(-max_d(:,i),'r','linewidth',2);
plot(-min_d(:,i),'r','linewidth',2);
title(['Mode water upper and lower limits :- PIES ' letra(i)])
hold off;
end;

print -depsc2 MW_min_max.eps
print -dpng MW_min_max.png

figure(2);
for i = 1:4
subtightplot(4,1,i,0.08,0.05,0.05);
plot(thick(:,i),'r','linewidth',2);
title(['Mode water layer thickness :- PIES ' letra(i)])
end;

save mw_thickness_vars.mat max_d min_d thick dep;

print -depsc2 MW_thickness.eps
print -dpng MW_thickness.png
