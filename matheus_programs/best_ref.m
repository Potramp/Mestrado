%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AFTER RUNNING full_run, THIS WILL CHOOSE THE BEST REFERENCE LEVEL 
% FOR THE SEASONAL CORRECTION.
%
%
%
%
%
%
%
% -CORTEZI, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

%%%% create a variable with all errors from the p_lvl folders (100, 120, 140, etc.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Make variable with all error from PIES and ISAS 16 deg depth...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foldernumber = [100:20:360];
isot16dif = NaN.*ones(length(foldernumber),72,4);
for i = 1:length(foldernumber)
	str = (['cd /Users/matheusvcortezi/Dropbox/Mestrado/matheus_programs/' num2str(foldernumber(i)) '/']);
	eval(str);
	load error;
	isot16dif(i,:,:) = error';	
	isot16depth_monthly(i,:,:) = monthly_depth;
end;

%%%%  SAVE THE NEW VARIABLE;
cd ..
save ~/Dropbox/Mestrado/matheus_programs/isot16dif.mat isot16dif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Find smalest diference among p levels...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%least_dif = find(nanmean(isot16dif,2) == min(nanmean(isot16dif,2)));

%%%% TAKE AWAY THE NaNs FOR THE CUMSUM
isot16dif(isnan(isot16dif)) = 0;

isot_summup = cumsum(isot16dif,2);
isot_summup = isot_summup(:,end,:);
isot_summup = squeeze(isot_summup);
save isot_summup.mat isot_summup
letra = ['ABCD']
figure
for i = 1:4;
        subtightplot(4,1,i,0.1,[0.1 0.05],[0.1 0.05]);
plot([100:20:360],squeeze(1-isot_summup(:,i)./mean(isot_summup(:,i))),'-*k','linewidth',3)
	title(['PIES ' letra(i)]);
	ylabel('Normalized Anomaly');
end;
	xlabel('seasonal correction reference level (dbar)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Plots...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
letra = ['ABCD']
in = 733864;
out = in + 83*30-24;
junk = in:34.25:out;
junk = junk(1:72);
for i = 1:4
	%subplot(4,1,i)
	subtightplot(4,1,i,0.1,[0.1 0.05],0.05);
	plot(junk,squeeze(-isot16depth_monthly(:,i,:))','k');
	title(['PIES ' letra(i)]);
	datetick('x','yyyy');
	axis 'tight';
	ylabel('depth (m)')
end;


print -dpng isot16.png


junk = cumsum(isot_summup,2);
find(junk(:,end) == max(junk(:,end)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('END OF LINE...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

