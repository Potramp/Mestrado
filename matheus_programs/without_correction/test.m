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

figure
for i = 1:4
	subtightplot(4,1,i,0.05,0.05,0.05);
	hold on;
	plot(timeaxis, -min_d(:,i), 'b');
	plot(timeaxis, -max_d(:,i), 'b');
	datetick('x','mmmyy');
	grid on;
end;

%%%% LOAD HYCOM DATA
%% HYCOM TRANSPORT
load ~/Dropbox/Mestrado/Hycom/bc_transp
%load /Users/'Matheus Cortezi'/Dropbox/Mestrado/Hycom/bc_transp
date(97) = [];
V(97) = [];

%-.-.-.-.-.-.-.-.-.-.-.-.- FIX DATES -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%% THE HYCOM TRANSPORT AND THE PIES TIME SERIES DO NOT START AND END AT
%%%% THE SAME POINT IN TIME. THIS HERE FIXES IT.

ltv = [date(1)  date(end); timeaxis(1) timeaxis(end)];

%%%% WE USE TIME STEPS OF 6 DAYS BECAUSE OF HYCOM
%%%% SO WE NEED TO FIND THE FIRST VALUE ON date THAT
%%%% IS A DATE AT LEAST SIX DAYS


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
imin_d = min_d;
imax_d = max_d;


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
        a = find(max_d(:,i) >= m1(i) +2.*s(i,1));
        b = find(min_d(:,i) >= m2(i) +2.*s(i,1));

        imax_d(a,i) = NaN;
        imin_d(b,i) = NaN;
        %max_d(a,i) = NaN;
        %min_d(b,i) = NaN;
        end;

figure
for i = 1:4
        subtightplot(4,1,i,0.05,0.05,0.05);
        hold on;
        plot(timeaxis, -min_d(:,i), 'ob');
        plot(timeaxis, -max_d(:,i), 'ob');
        plot(timeaxis, -imin_d(:,i) , 'r');
        plot(timeaxis, -imax_d(:,i), 'r');
	datetick('x','mmmyy');
        grid on;
end;



