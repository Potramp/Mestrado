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
close all; clear all;
more('off')

%%%% CHECK WHERE YOU ARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Current folder...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pwd

path = (['~/Dropbox/Mestrado/matheus_programs/without_correction/13_16_15']);


%%%% GET VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Load variables...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load w.mat;

%%%% time var. for plots
load ../timeaxis;

letra = ['ABCD'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. min_d, max_d, thick, Tt, Tb...')
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
                        Tt(j,i) = wT(v(1),j,i);
                        Tb(j,i) = wT(v(end),j,i);
                else;
                        min_d(j,i) = NaN;
                        max_d(j,i) = NaN;
                        Tt(j,i) = NaN;
                        Tb(j,i) = NaN;
		end;
        end;
end;

%%%% THICKNESS
thick = (max_d - min_d);
%thick(isnan(thick)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. Plots and figure')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% MIN DEPTH
figure(1)
for i = 1:4;
	subplot(4,1,i)
	plot(timeaxis,-min_d(:,i),'-'),datetick('x','yyyy')
	title(['Top :- PIES ' letra(i)])
	
	hold on;
	plot(timeaxis,-max_d(:,i),'-m'),datetick('x','yyyy')
	hold off;

	str = (['print -dpsc ' path '/top.eps']);
        eval(str);
end;

%%%% MAX DEPTH
figure(2)
%hold on;
for i = 1:4;
	subplot(4,1,i)
	plot(timeaxis,-max_d(:,i),'-m'),datetick('x','yyyy')
	title(['Bottom :- PIES ' letra(i)])
	str = (['print -dpsc ' path '/bottom.eps']);
        eval(str);

end;
%hold off;

%%%% LAYER THICKNESS
figure(3)
%figure(2)
for i = 1:4;
	subplot(4,1,i)
	plot(timeaxis,thick(:,i),'-'),datetick('x','yyyy')
	title(['Thickness :- PIES ' letra(i)])
	str = (['print -dpsc ' path '/thick.eps']);
        eval(str);
end;

%%%% Tt Tb
figure(4)
for i = 1:4;
	subplot(4,1,i)
	plot(timeaxis,Tt(:,i),'-'),datetick('x','yyyy')
	
	hold on;
	plot(timeaxis,Tb(:,i),'-m'),datetick('x','yyyy')
	hold off;

	title(['Tt and Tb :- PIES ' letra(i)])
	str = (['print -dpsc ' path '/TtTb.eps']);
        eval(str);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

