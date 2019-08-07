function [junk] = func_T_seasonal_model(my_p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This will calculate a seasonal model for T as a complement to 
% ds_find_T_from_tau.m.
%
% The idea is to follow Watts 2001 and make further seasonal corrections
% and see if that yields mode water formation.
%
% The theory as mathematics for this is in the paper. Read the appendix.
%
% 
%
% 
% CORTEZI 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd /Users/matheusvcortezi/Dropbox/Mestrado/matheus_programs/
junk = NaN;

%close all;
%clear all;
more('off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Loading variables...') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Tg_series
load CTD_and_Argo_p_levels
load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts/CleanUp_HydroData_ForGEM
%load /Users/apple/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts/CleanUp_HydroData_ForGEM


%%%% THIS PART MUST SELECT T VALUES UP TO THE SELECTED REFERENCE LEVEL FOR THE
%%%% TAU CORRECTION, ALWAYS CHECK THIS!

%%%% P level of reference. May be changed if necessary to test diferences.
%my_p = 200;
%my_p = 300;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['!!! YOUR P LEVEL OF REFERENCE IS ' num2str(my_p) ' db !!!'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indp = find(dep == my_p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Defining Tg...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tg = Tgseries(1:indp,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Defining T...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp=hydro(:,2);
pres=hydro(:,1);

%%%% The pressure levels of measures and interpolation are the same.
%%%% This selects each cas as a column, so Tg and T (in situ T) can be subtracted
%%%% as stated above.
junk = NaN.*ones(max(indxl-indxf),length(header));
for i = 1:length(indxf);
	junk(1:(indxl(i)-indxf(i)+1),i) = temp(indxf(i):indxl(i));
end;


%%%% Just a reminder that this will be used. 
%indp = find(dep == my_p);

T = junk(1:indp,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. Calculating T_prime...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_prime = NaN.*ones(indp,12);
for month = 1:12;
        indmon = find(header(:,5) == month);
        T_prime(:,month) = mean(T(:,indmon),2) - mean(Tg(:,indmon),2); 
end; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('5. Saving...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THIS WAS FOR 200
%save T_prime.mat T_prime

%%%% CHECKING THE EXISTENCE OF FILES AND FOLDERS WITH THE NAME OF my_p
folder = exist(['~/Dropbox/Meatrado/matheus_programs/' num2str(my_p)]);


%%%% CHEKING IF IT IS A FOLDER
if (folder ~= 7);
	%%%% MAKE THE FOLDER
	str = (['!mkdir ~/Dropbox/Meatrado/matheus_programs/' num2str(my_p)]);
	eval(str);
 
	%%%% SAVE THE VARIABLE
	str = (['save ~/Dropbox/Mestrado/matheus_programs/' num2str(my_p) '/T_prime.mat T_prime;']);
	eval(str);

elseif (folder == 7);	
	%%%% SAVE THE VARIABLE
	str = (['save ~/Dropbox/Mestrado/matheus_programs/' num2str(my_p) '/T_prime.mat T_prime;']);
	eval(str);
end;	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('6. Butterworth filter...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% FIXING MATRIX SIZE FOR ITERATION
indline = size(T_prime);
indline = indline(1);

for i = 1:indline %%%% number of lines
        str = (['rep_T_prime = repmat(T_prime,1,3);']);
        eval(str);
        [B,A] = butter(2,0.25);
        str = (['Y(i,:) = filtfilt(B,A,rep_T_prime(i,:));']);
        eval(str)
end;
T_tilde = Y(:,13:24);
%%%% And plot!
contourf(1:12,-dep(1:indline),Y(:,13:24))
%%%% And save!
str = (['save ~/Dropbox/Mestrado/matheus_programs/' num2str(my_p) '/T_tilde.mat T_tilde;']);
eval(str);
%%%% And print!
str = (['print -dpng ~/Dropbox/Mestrado/matheus_programs/' num2str(my_p) '/T_tilde.png;']);
eval(str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

