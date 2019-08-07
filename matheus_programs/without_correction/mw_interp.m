%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD mw_thickness_vars ON /Users/matheusvcortezi/Dropbox/Mestrado/matheus_programs/without_correction
% WE NEED TO INTERPOLATE SOME DATA TO FILL GAPS.
% THE DATA WITHOUT GAPS WILL HAVE IT'S CORRELATION CALCULATED WITH THE OUTPUT FROM
% transport.m 
%
%
% -CORTEZI,2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% clearing variables
clear all;
close all;
more('off');

%%%% LOADING MW LAYER VARIABLES
load mw_thickness_vars;
load timeaxis;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. finding NaNs...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inmax = isnan(max_d);
inmin = isnan(min_d);
imax_d = max_d;
imin_d = min_d;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Interp for dates with NaNs on trasport...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% SOME STATIONS HAVE DISCONTINUITIES ON DATA
%%%% AND THERE'S NEED TO ELIMINATE THOSE GAPS FROM INTERPOLATIONS
%%%% AND ALSO PREVENT EXTRAPOLATIONS

%%% DAYS OF GAPS
%dgap = [733850 NaN NaN 735883; 733899 734219 734691 736278; 733851 NaN NaN 735879; 733852 735484 735579 735879];
%tgap = NaN.*dgap;

%%% INDEXES OF GAPS, START AND END LIMITS;
%for i = 1:4
%	for j = 1:4;
%		if(~isnan(dgap(i,j)));
%			tgap(i,j) = find(timeaxis == dgap(i,j));
%		else;
%			tgap(i,j) = NaN;
%		end;
%	end;
%end;


%%%% INTERP MAX_D
%for i = 1:4
%	nmax = find(inmax(:,i) == 1);
%	numax = find(inmax(:,i) == 0);
%	numax = numax(end);
%	nmax = nmax(nmax<=numax);
%
%	%%%% GAPS ON STATIONS B AND D
%	if( i == 2 | i == 4);
%		nmax = nmax(nmax >= tgap(i,1));
%		nmax = nmax(nmax <= tgap(i,4));
%		nmax = nmax(nmax <= tgap(i,2) & nmax >= tgap(i,3));
%	else;
%		nmax = nmax(nmax >= tgap(i,1));
%		nmax = nmax(nmax <= tgap(i,4));
%	end;
%
%	yy = spline(timeaxis(~isnan(max_d(:,i))),max_d(~isnan(max_d(:,i)),i),timeaxis(nmax));
%	imax_d(nmax,i) = yy;
%end;

%%%% INTERP MIN_D
%for i = 1:4
%	nmin = find(inmin(:,i) == 1);
%	numin = find(inmin(:,i) == 0);
%	numin = numin(end);
%	nmin = nmin(nmin<=numin);
%
%	%%%% GAPS ON STATIONS B AND D
%	if( i == 2 | i == 4);
%		nmin = nmin(nmin >= tgap(i,1));
%		nmin = nmin(nmin <= tgap(i,4));
%		nmin = nmin(nmin <= tgap(i,2) & nmin >= tgap(i,3));
%	else;
%		nmin = nmin(nmin >= tgap(i,1));
%		nmin = nmin(nmin <= tgap(i,4));
%	end;
%
%	yy = spline(timeaxis(~isnan(min_d(:,i))),min_d(~isnan(min_d(:,i)),i),timeaxis(nmin));
%	imin_d(nmin,i) = yy;
%end;

%%%% without interpolation
figure(1)
for i = 1:4;
subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
hold on;
plot(timeaxis,max_d(:,i));
plot(timeaxis,min_d(:,i));
datetick('x','mmmyy');grid on;ylim([-200 600]);
end;
hold off;

print -depsc2 mw_max_min_d.eps

%i_thick = imax_d - imin_d;

%%%% interpolated
%figure(2)
%for i = 1:4
%subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
%hold on;
%plot(timeaxis,imax_d(:,i),'r','linewidth',2);
%plot(timeaxis,imin_d(:,i), 'b', 'linewidth', 2);
%plot(timeaxis,i_thick(:,i), 'k', 'linewidth', 2);
%datetick('x','mmmyy');grid on;ylim([-200 600]);
%end;
%hold off;

%print -depsc2 mw_interp_max_min_d.eps

%figure(3)
%for i = 1:4
%subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
%hold on;
%plot(timeaxis,i_thick(:,i), 'ok', 'linewidth', 2);
%plot(timeaxis,thick(:,i), 'r', 'linewidth', 2);
%datetick('x','mmmyy');grid on;ylim([-200 600]);
%end;
%hold off;

%print -depsc2 mw_thick.eps

%save interp_mw_thick.mat imax_d imin_d i_thick

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Resolution drop...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all;
%clear all;

%%%% WE SHOULD LOWER THE RESOLUTION TO MATCH THE ONE ON
%%%% THE HYCOM DATA. 

%%%% 2 METHODS MIGHT BE TESTED, FOR MEAN AND MEDIAN ALIKE.
%%%% GROUP DAYS -BEFORE- THE HYCOM DATE, AND -AFTER- IT.

%%%% LOAD HYCOM DATA 
%% HYCOM TRANSPORT
load ~/Dropbox/Mestrado/Hycom/bc_transp
%load /Users/'Matheus Cortezi'/Dropbox/Mestrado/Hycom/bc_transp
date(97) = [];
V(97) = [];
figure
plot(date,V);
hold on;
plot([date(1)  date(end)],[mean(V) mean(V)], 'm', 'linewidth', 2);
plot([date(1)  date(end)],[mean(V)+std(V) mean(V)+std(V)], 'm', 'linewidth', 2);
plot([date(1)  date(end)],[mean(V)-std(V) mean(V)-std(V)], 'm', 'linewidth', 2);
datetick('x','mmmyy'); grid on; 
xlabel('Time');
ylabel('Transport (Sv)')

mean(V)
std(V)

%-.-.-.-.-.-.-.-.-.-.-.-.- FIX DATES -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%% THE HYCOM TRANSPORT AND THE PIES TIME SERIES DO NOT START AND END AT
%%%% THE SAME POINT IN TIME. THIS HERE FIXES IT.

ltv = [date(1)  date(end); timeaxis(1) timeaxis(end)];

%%%% WE USE TIME STEPS OF 6 DAYS BECAUSE OF HYCOM
%%%% SO WE NEED TO FIND THE FIRST VALUE ON date THAT 
%%%% IS A DATE AT LEAST SIX DAYS BEFORE THE FIRST ONE ON timeaxis


%%%% start HAS THE FIRST VALUE WITHIN THE INTERVAL THAT COMPREHENDS PIES DATA
start = find(date >= timeaxis(1));
start = start(1);


%%%% WE WILL CALCULATE MEAN 6-DAY VALUES WITH DIFFERENT METHODS PRESENTED BELOW
%%%% FOR THAT WE NEED TO START AT DIFFERENT DAYS

%% ENDING POINT WHEN THE AVERAGE STARTS BEFORE THE DAY ON V
ending = find(date <= timeaxis(end));
ending = ending(end - 1);
%start_b(2) = 424;           %% +6 FALLS OUT OF RANGE

%% ENDING POINT WHEN THE AVERAGE STARTS AFTER THE DAY ON V
%start_a(2) = ending; 
%clear ending;

%% NOW WE NEED THE VALUES FOR THOSE DATES, TO LIMIT timeaxis
%tal = [date(start_b(1)) date(start_b(2)); date(start_a(1)) date(start_a(2))] 
tl = [start ending];

%%%% INDICES OF timeaxis DATES THAT ARE EQUAL TO date DATES
for i = start:ending;      
	a(i) = find(timeaxis == date(i));
end;
a = a(start:end);
a = a(2:end);
%-.-.-.-.-.-.-.-.-.-.-.-.- FIX DATES -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-


%-.-.-.-.-.-.-.-.-.-.-.-.- MEAN -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3.1. Mean before date...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% TURNING OUTLYERS INTO NANS
        m1 = nanmean(max_d);
        m2 = nanmean(min_d)
        for i = 1:4
                junk1  = max_d(:,i);
                junk2  = min_d(:,i);
                s(i,1) = std(junk1(~isnan(junk1)));
                s(i,2) = std(junk2(~isnan(junk2)));
        end;
        for i = 1:4
        c = find(max_d(:,i) <= m1(i) -2.*s(i,1));
        d = find(min_d(:,i) <= m2(i) -2.*s(i,1));

        max_d(c,i) = NaN;
        min_d(d,i) = NaN;
        end;
clear c d;


%%%% MATCH DATES
%Vmb = V(start_b(1):start_b(2));
Vmb = V(start+1:ending);
%timeb = timeaxis(start_b(1):start_b(2));
timeb = date(start+1:ending);


%%%% THE ACTUAL AVERAGE
%cont = 1;
for j = 1:4;
	for i = 1:length(a)
		if(i + 6 <= length(timeaxis))
			mbt_maxd(i,j) = nanmean(max_d(a(i)-5:a(i),j)); 
			mbt_mind(i,j) = nanmean(min_d(a(i)-5:a(i),j)); 
		end;
	end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3.2. Mean after date...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% MATCH DATES
%Vma = V(start+1:ending);
%Vma = V(start_a(1):start_a(2));
%timea = timeaxis(start_a(1):start_a(2));
%timea = date(start_a(1):start_a(2));
%timea = date(start+1:ending);

%%%% THE ACTUAL AVERAGE
%for j = 1:4;
%	for i = 1:length(a)
%		if(i + 6 <= length(timeaxis))
%			mat_maxd(i,j) = nanmean(max_d(a(i):5+a(i),j)); 
%			mat_mind(i,j) = nanmean(min_d(a(i):5+a(i),j)); 
%		end;
%	end;
%end;

%-.-.-.-.-.-.-.-.-.-.-.-.- MEAN -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

%-.-.-.-.-.-.-.-.-.-.-.-.- MEDIAAN -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('3.3. Median before date...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% THE ACTUAL MEDIAN
%for j = 1:4;
%	for i = 1:length(a)
%		if(i + 6 <= length(timeaxis))
%			mdbt_maxd(i,j) = nanmedian(max_d(a(i)-5:a(i),j)); 
%			mdbt_mind(i,j) = nanmedian(min_d(a(i)-5:a(i),j)); 
%		end;
%	end;
%end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('3.4. Median after date...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% THE ACTUAL MEDIAN
%for j = 1:4;
%	for i = 1:length(a)
%		if(i + 6 <= length(timeaxis))
%			mdat_maxd(i,j) = nanmedian(max_d(a(i):5+a(i),j)); 
%			mdat_mind(i,j) = nanmedian(min_d(a(i):5+a(i),j)); 
%		end;
%	end;
%end;

%-.-.-.-.-.-.-.-.-.-.-.-.- MEDIAAN -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%-.-.-.-.-.-.-.-.-.-.-.-.- LAYER THICKNESS -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. Mw thickness...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmbt = mbt_maxd - mbt_mind;
%tmat = mat_maxd - mat_mind;
%tmdbt = mdbt_maxd - mdbt_mind;
%tmdat = mdat_maxd - mdat_mind;

clear a;
clear i;
clear j;
clear m1 m2 s start thick tl ending V;
clear  V a dep i_thick i* n* str yy;

save mw_interp_vars.mat;

%-.-.-.-.-.-.-.-.-.-.-.-.- LAYER THICKNESS -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%-.-.-.-.-.-.-.-.-.-.-.-.- STATISTICS OF RESOLUTION DROPS -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('5. Statistics of resolution drops...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(12)
%list = ['tmbt ';'tmat ';'tmdbt';'tmdat'];
%list = ['tmbt'' ';'tmat'' ';'tmdbt''';'tmdat'''];
%for i = 1:4;
%	subtightplot(2,2,i,0.05,0.05,0.05)
%	hold on;
%	str = (['plot(1:4, nanmean(' list(i,:) '))']);
%	eval(str);
%	%str = (['plot(1:4, nanstd(' list(i,:) '), ''r'')']);
%	%eval(str);
%	hold off;
%end;

%-.-.-.-.-.-.-.-.-.-.-.-.- STATISTICS OF RESOLUTION DROPS -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%-.-.-.-.-.-.-.-.-.-.-.-.- INTERPOLATION -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('6. Interp low res...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% SOME STATIONS HAVE DISCONTINUITIES ON DATA
%%%% AND THERE'S NEED TO ELIMINATE THOSE GAPS FROM INTERPOLATIONS
%%%% AND ALSO PREVENT EXTRAPOLATIONS

%% CLEARING SOME VARIABLES
clear  V a dep i_thick i* n* str yy;



dgap = [733850 735883; 734691 736278; 733851 735878; 733870 735484];

%tnan = find(isnan(tmbt)); 
for i = 1:4
	%%%% COPY THE STATION TO BE INTERPOLATED
	junk = tmbt(:,i);
	str = (['pies' num2str(i) ' = junk(~isnan(junk));']);
	eval(str);
	%%%% DAYS OF NON-NAN DATA. USED ON INTERP AS THE t OF f(t)
	timejunk = timeb(~isnan(junk));

	%% FIND THE TIME WE REALLY WANT TO INTERPOLATE
	timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
	timelim =  timeb(timelim);

	%% INTERPOLATED VARIABLE. NO NANS INSIDE TIME INTERVAL
	str= (['i_tmbt' num2str(i) ' = interp1(timejunk, pies' num2str(i) ', timelim , ''linear'');']);
	eval(str);
	
end;


%-.-.-.-.-.-.-.-.-.-.-.-.- INTERPOLATION -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%-.-.-.-.-.-.-.-.-.-.-.-.- SINFIT AND FILTERING -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('7. Sinfit and filtering...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. SIN BEFORE BLACKMAN -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. 
%%%% sinfit
%% MW THICKNESS
for i = 1:4
	%% FIND THE TIME WE REALLY WANT TO USE
        timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
        t =  timeb(timelim);
	T = 365;
	str = (['coef = sinfitB(t, i_tmbt' num2str(i) ', T);']);
	eval(str);
	str = (['thicksin' num2str(i) '=coef(1)+coef(2)*t+coef(3)*sin((2*pi/T)*t)+coef(4)*cos((2*pi/T)*t);']);
	eval(str);
	str = (['thicksin' num2str(i) ' = thicksin' num2str(i) '(thicksin' num2str(i) ' ~=0);']);
	eval(str);
end;

% TRANSPORT
coef = sinfitB(timeb, Vmb, T);
%Vsin=coef(3)*sin((2*pi/T)*timeb)+coef(4)*cos((2*pi/T)*timeb);
Vsin=coef(1)+coef(2)*timeb+coef(3)*sin((2*pi/T)*timeb)+coef(4)*cos((2*pi/T)*timeb);
%Vsin = Vsin(Vsin ~=0);


%%% SUBTRACT SINFITS
% mw thickness
for i = 1:4;
	str = (['tmbt_dean' num2str(i) ' = - thicksin' num2str(i) ' + i_tmbt' num2str(i) ';']);
	eval(str);
end;

% transport
%V_dean = Vmb;
V_dean = Vmb - Vsin;

%%%% PLOT TO UNDERSTAND WHAT IS GOING ON;
figure(13); 
plot(Vsin, 'k'); 
hold on; 
plot(Vmb), 'b'; 
plot(V_dean,'r');
legend('Sinfit','Original','Minus Sinfit');
grid on;

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. SIN BEFORE BLACKMAN -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-. 

%%%% FILTERING WITH BLACKMAN
%% thickness
for j = 1:4


%%%% BLACKMAN WITH SINFIT SUBTRACTED
str = (['tmbt_black' num2str(j) ' = conv(tmbt_dean' num2str(j) ',blackman(31)/sum(blackman(31)),''same'');']);
eval(str);

%%%% BLACKMAN WITHOUT SINFIT SUBTRACTED
%str = (['tmbt_black' num2str(j) ' = conv(i_tmbt' num2str(j) ',blackman(31)/sum(blackman(31)),''same'');']);
%eval(str);
%str = (['tmbt_norm' num2str(j) ' = conv(ones(size(i_tmbt' num2str(j) ')),blackman(31)/sum(blackman(31)),''same'');']);
%eval(str);
%str = (['tmbt_black' num2str(j) ' = tmbt_black' num2str(j) './tmbt_norm' num2str(j) ';']);
%eval(str);

end;

%% TRANSPORT
% blackman after sinfit
V_black = conv(V_dean,blackman(31)/sum(blackman(31)), 'same' );


%-.-.-.-.-.-.-.-.-.-.-.-.- SINFIT AND FILTERING -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- CORRELATION -.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('8. Correlation...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% NOW WE HAVE MW_THICKNESS AND NEED TO CALCULATE THE CORRELATION WITH THE
%%%% BC TRANSPORT WE INTERPOLATED. 


for i = 1:4
        timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
	
	%%%% BLACKMAN not using SINFIT
	Vi = -V_black(timelim);
	%Vi =V_black(timelim);
	str = (['[c' num2str(i) ', lags' num2str(i) '] = xcov(Vi, tmbt_black' num2str(i) ', ''coef'');']);
	eval(str);
	str = (['[rho' num2str(i) ', p' num2str(i) '] = corrcoef(Vi, tmbt_black' num2str(i) ' );']);
	eval(str);
	
	% DO NOT USE THIS BIT
	%%%%%% BLACKMAN using SINFIT
	%%%Vi = V_dean(timelim);
	%%%str = (['[c' num2str(i) ', lags' num2str(i) '] = xcov(Vi, tmbt_dean' num2str(i) ', ''coef'');']);
	
end;


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- CORRELATION -.-.-.-.-.-.-.-.-.-.-.-.-.-

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.- PLOTS -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('9. Plots...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
letra = ['ABCD'];

figure(5)
for i = 1:4;
subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
if(i ==1)
	title('Mean before date')
end;
hold on;
plot(timeb,mbt_maxd(:,i));
plot(timeb,mbt_mind(:,i));
plot(timeb,tmbt(:,i),'k');
datetick('x','mmmyy');grid on;ylim([-200 600]);
grid on;
end;
hold off;

print -depsc2 1.eps

%figure(6)
%for i = 1:4;
%subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
%if(i ==1)
%title('Mean after date')
%end;
%hold on;
%plot(timea,mat_maxd(:,i));
%plot(timea,mat_mind(:,i));
%plot(timea,tmat(:,i),'k');
%datetick('x','mmmyy');grid on;ylim([-200 600]);
%end;
%hold off;
%print -depsc2 2.eps

%figure(7)
%for i = 1:4;
%subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
%if(i ==1)
%title('Median after date')
%end;
%hold on;
%plot(timeb,mdbt_maxd(:,i));
%plot(timeb,mdbt_mind(:,i));
%plot(timeb,tmdbt(:,i),'k');
%datetick('x','mmmyy');grid on;ylim([-200 600]);
%end;
%hold off;
%print -depsc2 3.eps

%figure(8)
%for i = 1:4;
%subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
%if(i ==1)
%title('Median after date')
%end;
%hold on;
%plot(timea,mdat_maxd(:,i));
%plot(timea,mdat_mind(:,i));
%plot(timea,tmdat(:,i),'k');
%datetick('x','mmmyy');grid on;ylim([-200 600]);
%end;
%hold off;
%print -depsc2 4.eps

figure(9)
for i = 1:4;
subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
hold on;
        timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
        timelim =  timeb(timelim);
	%t =  timeb(timelim);
str = (['plot(timelim,i_tmbt' num2str(i) ', ''linewidth'', 2);']);
eval(str)
str = (['plot([timelim(i) timelim(end)],[mean(i_tmbt' num2str(i) ') mean(i_tmbt' num2str(i) ')], ''m'', ''linewidth'', 2);']);
eval(str)
disp('MEAN AND STD')
str = ['mean(i_tmbt' num2str(i) ')'];
eval(str)
str = ['std(i_tmbt' num2str(i) ')'];
eval(str)
str = (['plot(timelim,thicksin' num2str(i) ',''r'', ''linewidth'', 2 );']);
eval(str)
datetick('x','mmyyyy');
grid on;
str = (['title(''PIES - ' letra(i) ''')']);
eval(str)
end;
hold off;
print -depsc2 sin_adjust.eps

figure(10)
for i = 1:4;
subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
        timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
        t =  timeb(timelim);
str = (['plot(t,tmbt_dean' num2str(i) ', ''linewidth'', 2);']);
eval(str)
datetick('x','mmyyyy');
grid on;
str = (['title(''PIES - ' letra(i) ''')']);
eval(str)
end;

figure(11)
for i = 1:4;
subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
t =  timeb(timelim);
str = (['plot(t,tmbt_black' num2str(i) ', ''linewidth'', 2);']);
eval(str)
%hold on;
%plot(timeb,V_black, 'k', 'linewidth', 2);
datetick('x','mmyyyy');
grid on;
ylabel('Layer thickness anomaly (m)')
str = (['title(''PIES - ' letra(i) ''')']);
eval(str)
end;
print -depsc2 blackman_filtered.eps

figure(12)
for i = 1:4;
subtightplot(2,2,i,0.1,[0.1 0.05],[0.1 0.05]);
str = (['plot(lags' num2str(i) ',c' num2str(i) ', ''linewidth'', 2);']);
eval(str)

%dim = [.2*i .5 .3 .3];
%str = (['p = ' num2str(p1)]);
%annotation('textbox',dim,'String',str,'FitBoxToText','on');

ylabel('Correlation');
xlim([-90 90])
grid on;
str = (['title(''PIES - ' letra(i) ''')']);
eval(str)
end;
xlabel('Lag')
print -depsc2 correlation.eps

figure(14)
plot(timeb,V_dean);
hold on;
plot(timeb,V_black, 'k', 'linewidth', 2);
plot([timeb(1)  timeb(end)],[mean(V_black) mean(V_black)], 'm', 'linewidth', 2);
grid on; datetick('x','yyyy');
xlabel('Time');
ylabel('Transport (Sv)')
hold off;
print -depsc -r150 BCtransp.eps

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.- PLOTS -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-


save covs.mat c* lags* V_black *sin dgap timeb T thicksin* date tmbt_black*
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
