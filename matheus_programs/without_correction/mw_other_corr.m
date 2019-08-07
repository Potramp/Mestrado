%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.TEST CORRELATION OF:
% -MW MEAN DEPTH AND TRANSPORT
% 	-with blackman
%	-without blackman
% -VSQUARED
% 2.COMPARE SEASONAL CYCLES OF: 
% -MAX_D, MIN_D, THICK, V.
%
%
%
% -CORTEZI 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
more('off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Loading variables...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% mw variables
load mw_interp_vars
load itmbt;

%%%% define temporal limits for each station. A B C D, in order. Start and end.
dgap = [733850 735883; 734691 736278; 733851 735878; 733870 735484];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Interpolations to fill gaps...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% THIS IS A COPY OF THE SAME SNIPET ON mw_interp, ONLY THAT IT IS MODIFIED TO
%%%% INTERPOLATE MAX AND MIN DEPTHS OF THE MODE WATER COLUMN
for i = 1:4
        %%%% COPY THE STATION TO BE INTERPOLATED
        temp1 = mbt_maxd(:,i);
        temp2 = mbt_mind(:,i);
        str = (['max' num2str(i) ' = temp1(~isnan(temp1) & ~isnan(temp2));']);
        eval(str);
        str = (['min' num2str(i) '  = temp2(~isnan(temp1) & ~isnan(temp2));']);
        eval(str);


        %%%% DAYS OF NON-NAN DATA. USED ON INTERP AS THE t OF f(t)
        timejunk = timeb(~isnan(temp1));

        %% FIND THE TIME WE REALLY WANT TO INTERPOLATE
        timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
        timelim =  timeb(timelim);

        %% INTERPOLATED VARIABLE. NO NANS INSIDE TIME INTERVAL
        str= (['imin' num2str(i) ' = interp1(timejunk, min' num2str(i) ', timelim , ''linear'');']);
        eval(str);
        str= (['imax' num2str(i) ' = interp1(timejunk, max' num2str(i) ', timelim , ''linear'');']);
        eval(str);


end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Mean depth of mw...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:4
	str = (['avd' num2str(i) ' = (imin' num2str(i) ' + imax' num2str(i) ')./2;']);
	eval(str);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. Adjusted seasonal cycles...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:4
        %% FIND THE TIME WE REALLY WANT TO USE
        timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
        t =  timeb(timelim);
        T = 365;

	%% THICKNESS
        %% FIND THE TIME WE REALLY WANT TO USE
        str = (['coef = sinfitB(t, i_tmbt' num2str(i) ', T);']);
        eval(str);
        str = (['thicksin' num2str(i) '=coef(1)+coef(2)*t+coef(3)*sin((2*pi/T)*t)+coef(4)*cos((2*pi/T)*t);']);
        eval(str);
        str = (['thicksin' num2str(i) ' = thicksin' num2str(i) '(thicksin' num2str(i) ' ~=0);']);
        eval(str);

	%% MIN DEPTH
        str = (['coef = sinfitB(t, imin' num2str(i) ', T);']);
        eval(str);
        str = (['minsin' num2str(i) '=coef(1)+coef(2)*t+coef(3)*sin((2*pi/T)*t)+coef(4)*cos((2*pi/T)*t);']);
        eval(str);
        str = (['minsin' num2str(i) ' = minsin' num2str(i) '(minsin' num2str(i) ' ~=0);']);
        eval(str);

	%% MAX DEPTH
        str = (['coef = sinfitB(t, imax' num2str(i) ', T);']);
        eval(str);
        str = (['maxsin' num2str(i) '=coef(1)+coef(2)*t+coef(3)*sin((2*pi/T)*t)+coef(4)*cos((2*pi/T)*t);']);
        eval(str);
        str = (['maxsin' num2str(i) ' = maxsin' num2str(i) '(maxsin' num2str(i) ' ~=0);']);
        eval(str);
	

	%% AVG DEPTH
        str = (['coef = sinfitB(t, avd' num2str(i) ', T);']);
        eval(str);
        str = (['avdsin' num2str(i) '=coef(1)+coef(2)*t+coef(3)*sin((2*pi/T)*t)+coef(4)*cos((2*pi/T)*t);']);
        eval(str);
        str = (['avdsin' num2str(i) ' = avdsin' num2str(i) '(avdsin' num2str(i) ' ~=0);']);
        eval(str);
	
end;

% TRANSPORT
coef = sinfitB(timeb, Vmb, T);
Vsin=coef(1)+coef(2)*timeb+coef(3)*sin((2*pi/T)*timeb)+coef(4)*cos((2*pi/T)*timeb);

%%%% SUBTRACTING ADJUSTED SINUSOID
for i = 1:4
	str = (['min_dean' num2str(i) ' = imin' num2str(i) ' - minsin' num2str(i) ';']);
	eval(str);
	str = (['max_dean' num2str(i) ' = imax' num2str(i) ' - maxsin' num2str(i) ';']);
	eval(str);
end

Vdean = Vmb - Vsin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('5. Blackman filter...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:4
	
%%%% BLACKMAN WITHOUT SINFIT SUBTRACTED
str = (['min_black' num2str(j) ' = conv(imin' num2str(j) ',blackman(65)/sum(blackman(65)),''same'');']);
eval(str);
str = (['min_norm' num2str(j) ' = conv(ones(size(imin' num2str(j) ')),blackman(65)/sum(blackman(65)),''same'');']);
eval(str);
str = (['min_black' num2str(j) ' = min_black' num2str(j) './min_norm' num2str(j) ';']);
eval(str);

str = (['max_black' num2str(j) ' = conv(imax' num2str(j) ',blackman(65)/sum(blackman(65)),''same'');']);
eval(str);
str = (['max_norm' num2str(j) ' = conv(ones(size(imax' num2str(j) ')),blackman(65)/sum(blackman(65)),''same'');']);
eval(str);
str = (['max_black' num2str(j) ' = max_black' num2str(j) './max_norm' num2str(j) ';']);
eval(str);

%%%% BLACKMAN FILTER WITH SINFIT SUBTRACTED
str = (['min_black_sin' num2str(j) ' = conv(min_dean' num2str(j) ',blackman(65)/sum(blackman(65)),''same'');']);
eval(str);
str = (['max_black_sin' num2str(j) ' = conv(max_dean' num2str(j) ',blackman(65)/sum(blackman(65)),''same'');']);
eval(str);

end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('6. Compare seasonal cycles...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SOME PLOTS FOR VISUAL COMPARISON
figure(1)
for i = 1:4;
	subtightplot(4,1,i,0.1,[0.1 0.05],0.05);
		timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
		t =  timeb(timelim);
		%t =  1:12;
	%% THICKNESS
	str = (['plot(t,(thicksin' num2str(i) '- mean(thicksin' num2str(i) '))./mean(thicksin' num2str(i) '), ''g'')']);
	eval(str);
	hold on;
	
	%% MIN DEPTH
	str = (['plot(t,(minsin' num2str(i) '- mean(minsin' num2str(i) '))./mean(minsin' num2str(i) '))']);
	eval(str);
	
	%% MAX DEPTH
	str = (['plot(t,(maxsin' num2str(i) '- mean(maxsin' num2str(i) '))./mean(maxsin' num2str(i) '), ''m'')']);
	eval(str);

	%% AVERAGE DEPTH
	str = (['plot(t,(avdsin' num2str(i) '- mean(avdsin' num2str(i) '))./mean(avdsin' num2str(i) '), ''k'')']);
	eval(str);
	

	%% TRANSPORT
	str = (['plot(timeb,-(Vsin-mean(Vsin))./mean(Vsin),''r'')']);
	eval(str);
	
	datetick('x','mmyyyy');
	xlabel('Date')
	ylabel('Normalized signal')
	grid on;
	legend('thickness','min','max','avg', 'cb');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('6. Correlation of average depth and transport seasonal cycles...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:4
        timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));

        %%%% BLACKMAN not using SINFIT
        Vi = Vsin(timelim);
        str = (['[c' num2str(i) ', lags' num2str(i) '] = xcov(-Vi, avdsin' num2str(i) ', ''coef'');']);
        eval(str);
        str = (['[d' num2str(i) ', dlags' num2str(i) '] = xcov(-Vi, thicksin' num2str(i) ', ''coef'');']);
        eval(str);

end;

figure(2)
for i = 1:4;
        subtightplot(4,1,i,0.1,[0.1 0.05],0.05);
	str = (['plot(lags' num2str(i) ',c' num2str(i) ', ''linewidth'', 2);']);
	eval(str);
	hold on;
	str = (['plot(dlags' num2str(i) ',d' num2str(i) ',''r'', ''linewidth'', 2);']);
	eval(str)
	xlabel('Lag')
	ylabel('Cross-covariance')
	legend('Average','Thickness');
end;
