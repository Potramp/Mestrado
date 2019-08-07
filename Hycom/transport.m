%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates transport for the BC near the PIES SAM array
%
%
%
%
%
% - CORTEZI, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
more('off')

%%%% LOAD DATA
 
load ~/Dropbox/Mestrado/Hycom/hycom_34S.mat;

%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Ploting V to find BC...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% DATE VECTOR
date = datenum(yy,mm,dd);

ilo = find(lon >= -52);
%lon = lon(ilo)
%vv = vv(:,ilo,:);

%%%% FINDING BC CORE
bc = [1 7 12];
for i = [1 2 3];
	subplot(1,3,i);

	%%%% V OVER TIME CONTOURS
	contourf(lon,date,squeeze(vv((bc(i)),:,:))',15);
	datetick('y','yyyy')
	hold on;

	%%%% V = 0 CONTOUR
	contour(lon,date,squeeze(vv(bc(i),:,:))',[0 0], 'linecolor', 'k', 'linewidth', 3);
	
	%%%% SUBPLOT TITLE
	str = (['title(''' num2str(dep(bc(i))) ' m '')']);
	eval(str);
	hold off;
	grid on;
end;
%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.
%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Averaging V...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% DEFINING VARIABLES TO SELECT DATA MONTHLY
ylim = [yy(1),yy(end)];
ylist = [ylim(1):ylim(2)];
%mlist = [1:12];


%%%% MONTHLY AVERAGE
cont = 1;
%% GO OVER YEARS
for i = 1:length(ylist);
	
	indy = find(yy == ylist(i));
	%% GO OVER MONTHS
	
	for j = 1:12; 
		indm = find(mm == j);
		
		%%%% RIGHT YEAR AND RIGHT MONTH
		indval = intersect(indy,indm);
		vmm(:,:,cont) = mean(vv(:,:,indval),3);
		cont = cont+1;
	end;
end;

%%%% fixing date vector
ymm = repmat(ylist,12);
ymm = (ymm(1:84));
%ymm = num2str(ymm(1:84));
mmm = repmat([1:12],7);
mmm = (mmm(1,:)); 
%mmm = num2str(mmm(1,:)); 
dmm = ones(1,length(ymm));
datemm = datenum(ymm,mmm,dmm);    
 
%%%% PLOTING MONTHLY AVERAGED V

figure
bc = [1 4 7];
for i = [1 2 3];
        subplot(1,3,i);

        %%%% V OVER TIME CONTOURS
        contourf(lon,datemm,squeeze(vmm((bc(i)),:,:))',15);
        datetick('y','yyyy')
        hold on;

        %%%% V = 0 CONTOUR
        contourf(lon,datemm,squeeze(vmm(bc(i),:,:))',[0 0], 'linecolor', 'k', 'linewidth', 3);

        %%%% SUBPLOT TITLE
        str = (['title(''' num2str(dep(bc(i))) ' m '')']);
        eval(str);
        hold off;
	grid on;
end;

clear ymm mmm dmm datemm;

%%%% MEAN V
v = mean(vmm,3);
	%plot
	figure
	
	%%%%% dep lim = 1000 m
	%contourf(lon,-dep(1:26),v(1:26,:));
	%%contourf(lon,-dep(1:26),v(1:26,:),15);
	%hold on;
	%contour(lon,-dep(1:26),v(1:26,:),[0 0], 'linecolor', 'k', 'linewidth', 3);
	%contour(lon,-dep(1:26),v(1:26,:),[-0.06 -0.06], 'linecolor', 'k', 'linewidth', 3);
	%colormap 'gray'
	%colorbar;
	%grid on;

	%%%% dep lim = 1000 m
	contourf(lon,-dep,v);
	%contourf(lon,-dep,v,15);
	hold on;
	contour(lon,-dep,v,[0 0], 'linecolor', 'k', 'linewidth', 3);
	contour(lon,-dep,v,[-0.06 -0.06], 'linecolor', 'k', 'linewidth', 3);
	
	%%%% PLOT SITES A TO D

	plot([-51.5 -49.5 -47.5 -44.5],[-dep(end) -dep(end) -dep(end) -dep(end)], 'gs', ...
		'LineWidth',2,...
    		'MarkerSize',10,...
    		'MarkerEdgeColor','b',...
    		'MarkerFaceColor',[0.5,0.5,0.5])
	colormap 'gray'
	colorbar;
	grid on;
	


%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('3. Check BC...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CHECK BC

%VALID INDEXES
%indval = find(lon <= -51.2 & lon >= );
indval = find(lon <= -51.2);

%%%% LOWER LIMIT FOR BC. THIS COMES FROM THE IMAGE FOR MEAN VELOCITIES
bcll = find(dep == 1500);

%%%% V LIMITED TO THE LENGTH AND DEPTH OF BC
bcv = vv(:,indval,:);
bcv = bcv(1:bcll,:,:);

%%%% MEAN V ON BC REGION
bcvm = mean(bcv,3);

%%%% PLOT MEAN V ON BC REGION
figure
contourf(lon(indval),-dep(1:bcll),bcvm); colormap 'gray'; 
	grid on;

clear indval bcv bcvm bcll;
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.
%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('4. Vertically integrated v at all longitudes...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% find lon of the first positive value for v_transp

%%%% finding nans
indnan = isnan(vv);

%%%% forced depth limit
depl = 21; % 1200 m
%depl = 17; % 800 m

%%%% CALCULATING TRANSPORT
for i = 1:length(date);
	for j = 1:length(lon);
		
		%%%% LAST VALID Z TO CALCULATE TRANSPORT
		%ll = find(indnan(:,j,i) == 1);
		%if(length(ll) == 0);
		%	depl = dep(end);
		%else;
		%	depl = dep(ll(1));
		%end;
		
		vtransp(j,i) = nansum(vv(1:depl,j,i).*dep(1:depl)');		
		%vtransp(j,i) = vtransp(j,i);
	end;
end;

%%%% PLOT HOVMOLLER OF TRANSPORT
%contourf(lon,date,vtransp'); datetick('y','yyyy'); colorbar; colormap 'gray';
%hold on; 
%contour(lon,date,vtransp',[0 0], 'linecolor', 'k', 'linewidth', 3); 
%hold off;

%%%% MEAN TRANSPORT
%vtm = nanmean(vtransp,2);
%figure
%plot(lon,vtm);


%%%% THE UNITS OF TRANSPORT ARE ALL MESSED UP. TO GET SVERDRUPS WE DO:
	%% width of transport section
	dx= diff(lon(1:end)*111000*cosd(34.5));

	%% equal the points on both variables
	half_vtransp = (vtransp(1:end-1,:) + vtransp(2:end,:))./2;

	%% fix units
	dx = repmat(dx,1,424);
	sv_vt = half_vtransp.*dx/1e6;

	%% clearing variables
	clear  half_vtransp;

%%%% PLOT HOVMOLLER OF TRANSPORT WITH THE RIGHT UNITS
%%fix lon size
sv_lon = (lon(1:end-1)+lon(2:end))./2;
contourf(sv_lon,date,sv_vt'); datetick('y','yyyy'); colorbar; colormap 'gray';
hold on; 
contour(sv_lon,date,sv_vt',[0 0], 'linecolor', 'k', 'linewidth', 3); 
	grid on;
hold off;

%%%% MEAN TRANSPORT ON SVERDRUPS
vtm = nanmean(sv_vt,2);
figure
plot(sv_lon,vtm);

%%%% SAVE RESULTS AND CLEAR VARIABLES
save sv_transport.mat sv_vt sv_lon dx date;
clear sv_vt sv_lon dx;


%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.
%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('5. Calculate transport for BC...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% LON AND LAT CAME AS SINGLE, MAKING IT'S TYPE DOUBLE
lon = double(lon);
lat = double(lat);

%%%% BC TRANSPORT
%% MAXIMUM WEST LIMIT
uwl = ilo(1);  % -52
%% MINIMUM WEST LIMIT
lwl = 32;  % -51.52

%% MINIMUM EAST LIMIT
lel = 40; % -50.5601
%% MAXIMUM EAST LIMIT
uel = 70; % -48.4800

clear dx;
for i = 1:length(date)
	%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
	%%%% WESTERN LIMIT
	negs = find(vtransp(:,i) <=0);
	el = negs(end);
	if(negs(1) <= uwl )
		wl = uwl;
	elseif(negs(1) >= lwl)
		wl = lwl;
	else;
		wl = negs(1);
	end;	
	%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
	
	%find were wl skips numbers
	for k = 1:length(negs)-1;
		if(negs(k+1) ~= negs(k)+1)
                       el = negs(k);
		end;
	end;	

	%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
	%%%% EASTERN LIMIT
	if(el >= uel )
                el = uel;
        elseif(el <= lwl)
                el = lel;
        end;	
	%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
	valids = (lon<= lon(el) & lon >= lon(wl));
	bcti(:,i) = vtransp(:,i).*valids;
	
	
	%%%% CB LENGTH. dx HAS DISTANCE BETWEEN BC LIMITS.
	dx(i) = pos2dist(lat,lon(wl),lat,lon(el),2);
	
end;

%%%% BCTI HAS VERTICAL INTEGRATED VELOCITIES. I NEED TRANSPORT. CONVERTING UNITS.
%% SUM O INTEGRATED V
bcti = sum(bcti,1);
bcti = bcti.*dx;

%%%% TRANSPORT IN SVERDRUPS
sv_bcti = bcti./1000000;

%%%% PLOT TRANSPORT TIME SERIES
figure;
plot(date,sv_bcti,'k','linewidth',2);
title('Brazil Current Transport');
ylabel('Transport (Sv)');
xlabel('t');
datetick('x','mmmyy');
grid on;

V = sv_bcti;
save bc_transp.mat V date;

print -depsc2 bctransp_sv.eps

%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.
%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.




