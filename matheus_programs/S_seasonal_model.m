%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This will calculate a seasonal model for S as a complement to 
% ds_find_S_from_tau.m.
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

close all;
clear all;
more('off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Loading variables...') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Sg_series
load CTD_and_Argo_p_levels
load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/Gem_data_and_scripts/CleanUp_HydroData_ForGEM


%%%% THIS PART MUST SELECT S VALUES UP TO THE SELECTED REFERENCE LEVEL FOR THE
%%%% TAU CORRECTION, ALWAYS CHECK THIS!

%%%% P level of reference. May be changed if necessary to test diferences.
my_p = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['!!! YOUR P LEVEL OF REFERENCE IS ' num2str(my_p) ' db !!!'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indp = find(dep == my_p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Defining Sg...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sg = Sgseries(1:indp,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Defining S...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sal=hydro(:,3);
pres=hydro(:,1);

%%%% The pressure levels of measures and interpolation are the same.
%%%% This selects each cas as a column, so Sg and S (in situ S) can be subtracted
%%%% as stated above.
junk = NaN.*ones(max(indxl-indxf),length(header));
for i = 1:length(indxf);
	junk(1:(indxl(i)-indxf(i)+1),i) = sal(indxf(i):indxl(i));
end;


%%%% Just a reminder that this will be used. 
%indp = find(dep == my_p);

S = junk(1:indp,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. Calculating S_prime...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_prime = NaN.*ones(11,12);
for month = 1:12;
        indmon = find(header(:,5) == month);
        S_prime(:,month) = mean(S(:,indmon),2) - mean(Sg(:,indmon),2); 
end; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('5. Saving...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save S_prime.mat S_prime


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('6. Butterworth filter...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:11 %%%% number of lines
        str = (['rep_S_prime = repmat(S_prime,1,3);']);
        eval(str);
        [B,A] = butter(2,0.25);
        str = (['Y(i,:) = filtfilt(B,A,rep_S_prime(i,:));']);
        eval(str)
end;
S_tilde = Y(:,13:24);
%%%% And plot!
contourf(1:12,-dep(1:11),Y(:,13:24))
%%%% And save!
save S_tilde.mat S_tilde;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

