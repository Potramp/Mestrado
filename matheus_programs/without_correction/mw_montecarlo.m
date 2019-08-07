%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONTECARLO SIMULATIONS FOR P-VALUE AND CORRELATION
%
%
%
%
%
%
%
%
%
% -CORTEZI, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
more('off')

%% loading data
load covs;

letra = ['ABCD'];

%%%% MONTECARLO SIMULATIONS TO GET P-VALUE

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%% EXPLORATORY STATS
for i = 1:4
        str = (['m_thick' num2str(i) ' = mean(tmbt_black' num2str(i) ');']);
        eval(str);
        str = (['s_thick' num2str(i) ' = std(tmbt_black' num2str(i) ');']);
        eval(str);
end;

m_V = mean(-V_black);
s_V = std(-V_black);
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%% STARTING VARIABLES   
for i = 1:4
        timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
        t =  timeb(timelim);
        t = length(t);
	str = (['density' num2str(i) ' = 0.*ones(1,2*length(t)-1);']);
	eval(str);        
end;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

%% CLEARING LAG VARIABLES
clear lags*;

%% MONTECARLO
%for i = 1:4
%	timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
%        t =  timeb(timelim);
%
%	for j = 1:10000
%                %t = length(t);
%                str = (['r = randn(size(tmbt_black' num2str(i) '))*std(tmbt_black' num2str(i) ');']);
%                eval(str)
%		rV = randn(size(V_black))*std(V_black);
%
%                %% FIT A SEASONAL SINUSOID
%                T = 365;
%                str = (['coef = sinfitB(t, r, T);']);
%                eval(str);
%                str = (['rsin=coef(1)+coef(2)*t+coef(3)*sin((2*pi/T)*t)+coef(4)*cos((2*pi/T)*t);']);
%                eval(str);
%                str = (['rsin = rsin(rsin ~=0);']);
%                eval(str);
%
%                coef = sinfitB(timeb, rV, T);
%                rVsin=coef(1)+coef(2)*timeb+coef(3)*sin((2*pi/T)*timeb)+coef(4)*cos((2*pi/T)*timeb);
%
%
%                %%%% REMOVE SEASONAL SIGNAL
%                r = r-rsin;
%                rV = rV - rVsin;
%
%                %%%% BLACKMAN FILTER
%                r = conv(r,blackman(31)/sum(blackman(31)),'same');
%                rV = conv(rV,blackman(31)/sum(blackman(31)),'same');
%
%		%%%% COV
%                str = (['[d' num2str(i) ',lags' num2str(i) '] = xcov(rV(timelim),r,''coef'');']);
%                eval(str)
%
%                str = (['density' num2str(i) ' = density' num2str(i) ' + (abs(d' num2str(i) ') >= abs(c' num2str(i) '));']); 
%                eval(str)
%        end;
%        str = (['density' num2str(i) ' = density' num2str(i) './10000;']);
%        eval(str);
%        
%end;

%%%% OLGA'S SUGGESTION

clear r*;

for i = 1:4
        for j = 1:10000
                timelim = find(timeb >= dgap(i,1) & timeb <= dgap(i,2));
                t =  timeb(timelim);
                %t = length(t);
                str = (['r = randn(size(tmbt_black' num2str(i) '))*std(tmbt_black' num2str(i) ');']);
                eval(str)
		rV = randn(size(V_black))*std(V_black);

		
		%% FIT A SEASONAL SINUSOID
	        T = 365;
	        str = (['coef = sinfitB(t, r, T);']);
	        eval(str);
	        str = (['rsin=coef(1)+coef(2)*t+coef(3)*sin((2*pi/T)*t)+coef(4)*cos((2*pi/T)*t);']);
	        eval(str);
	        str = (['rsin = rsin(rsin ~=0);']);
	        eval(str);
	
		coef = sinfitB(timeb, rV, T);
		rVsin=coef(1)+coef(2)*timeb+coef(3)*sin((2*pi/T)*timeb)+coef(4)*cos((2*pi/T)*timeb);


		%%%% REMOVE SEASONAL SIGNAL
		r = r-rsin;
		rV = rV - rVsin;

		%%%% BLACKMAN FILTER
		r = conv(r,blackman(31)/sum(blackman(31)),'same');
		rV = conv(rV,blackman(31)/sum(blackman(31)),'same');

                %%%% COV
                str = (['[d' num2str(i) ',lags' num2str(i) '] = xcov(rV(timelim),r,''coef'');']);
                eval(str)
		str = (['junk = abs(d' num2str(i) ');']);
                eval(str)
                %str = (['maxcov' num2str(i) '(j, 1, :) = abs(d' num2str(i) ') == max(abs(d' num2str(i) '));']);
                %eval(str)

		str = (['[pks locs] = findpeaks(-d' num2str(i) '(length(t):end));']);
                eval(str)
		
		%%%% GET THE POSITIVE VALUES
		pospks = pks(pks>0);
		poslocs = locs(pks>0);
		pk = pospks(1);
		loc = poslocs(1);
		
                str = (['maxcov' num2str(i) '(j, 2) = pk;']);
                eval(str)
                str = (['maxcov' num2str(i) '(j, 1) = loc;']);
                eval(str)
                %str = (['maxcov' num2str(i) '(j, 2) = max(junk);']);
                %eval(str)
                %str = (['maxcov' num2str(i) '(j, 1) = find(junk == max(junk));']);
                %eval(str)
                
		%%%% PVALUE
		str = (['density' num2str(i) ' = density' num2str(i) ' + (abs(d' num2str(i) ') >= abs(c' num2str(i) '));']); 
                eval(str)
        end;
        str = (['density' num2str(i) ' = density' num2str(i) './10000;']);
        eval(str);
end;




figure
for i = 1:4;
subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
str = (['plot(maxcov' num2str(i) '(:,1))']);
eval(str)
hold on;
str = (['title(''PIES - ' letra(i) ''')']);
eval(str)
end;

figure
for i = 1:4;
subtightplot(2,2,i,0.1,[0.1 0.05],0.05);
str  = (['plot(maxcov' num2str(i) '(:,2),''.r'')']);
eval(str)
str = (['title(''PIES - ' letra(i) ''')']);
eval(str)
end;
print -dpng -r150 ./blackman/p_cov.png

figure
for i = 1:4;
subtightplot(2,2,i,0.13,[0.1 0.05],[0.1 0.05]);
str = (['plot(6*lags' num2str(i) ', density' num2str(i) ',''r'',''linewidth'',2)']);
eval(str);
hold on;
str = (['plot(6*lags' num2str(i) ',c' num2str(i) ', ''linewidth'', 2);']);
eval(str)

%dim = [.2*i .5 .3 .3];
%str = (['p = ' num2str(p1)]);
%annotation('textbox',dim,'String',str,'FitBoxToText','on');

ylabel('Covariance');
xlim([0 360])
grid on;
str = (['title(''PIES - ' letra(i) ''')']);
eval(str)
xlabel('Lag')
print -dpng correlation.png
xlabel('Lag (days)')
%ylabel('P-value');
end;
print -dpng -r150 ./blackman/p_cov.png

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

