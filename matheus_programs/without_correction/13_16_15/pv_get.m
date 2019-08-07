%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I have to find the right pv to call mode water mode water.
% let's do it.
%
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
disp('3. Plot PV...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spv = size(PV);

for i = 1:spv(2)
	for j = 1:spv(3) 
		plot(PV(:,i,j),-mean_z(:,1,1));
		hold on;
	end;
end;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. Mean pv prof...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sT = size(wT);

for i = 1:spv(1)
	for j=1:spv(3)
		row = PV(i,:,j);
		notnan = row(~isnan(row));
		junk(i,j) = mean(notnan); 
	end;
end;

figure
plot(junk,squeeze(-mean_z(:,1,1)))
