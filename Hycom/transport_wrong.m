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
	
	%%%% dep lim = 1000 m
	contourf(lon,-dep(1:26),v(1:26,:));
	%contourf(lon,-dep(1:26),v(1:26,:),15);
	hold on;
	contour(lon,-dep(1:26),v(1:26,:),[0 0], 'linecolor', 'k', 'linewidth', 3);
	contour(lon,-dep(1:26),v(1:26,:),[-0.06 -0.06], 'linecolor', 'k', 'linewidth', 3);
	colormap 'gray'
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
bcll = find(dep == 500);

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
disp('4. Calculate transport at all longitudes...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% find lon of the first positive value for v_transp

%%%% finding nans
indnan = isnan(vv);


%%%% CALCULATING TRANSPORT
for i = 1:length(date);
	for j = 1:length(lon);
		ll = find(indnan(:,j,i) == 1);

		%%%% LAST VALID Z TO CALCULATE TRANSPORT
		if(length(ll) == 0);
			depl = dep(end);
		else;
			depl = dep(ll(1));
		end;
		vtransp(j,i) = nansum(vv(:,j,i));	
		vtransp(j,i) = vtransp(j,i).*depl;		
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

%%%% BC TRANSPORT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THIS DOESN'T WORK BUT IT'S TOO BEATIFUL FOR ME TO ERASE IT. AND IT TOOK A LOT OF WORK.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i = 1:length(date)
%	
%	%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%	%% FIND BC ZONAL LIMITS
%
%	%% WESTERN LIMIT ON CONTINENTAL SHELF BREAK (ACTUALLY A LITTLE FURTHER, TO GET AT LEAST 2000 m DEPTH)
%	%wl = find(lon >= -51.6799);
%	wl = find(vv(1,ilo,i) <= 0);
%	
%	%% THIS WILL TELL IF THE CONDITION OF IF WAS EVER REACHED;
%	cont = 0;
%	for k = 1:length(wl)-1
%		if(wl(k+1) ~= wl(k)+1)
%			el = wl(k);
%			cont = cont+1;
%			possible(cont) = k;
%		end;
%		
%	end;
%	
%
%	%%%% SOMETIMES THE CURRENT IS TO FAR EAST, SO OTHER MINOR NEGATIVE VELOCITIES CAUSE PROBLEMS, THIS SOLVES IT
%	%%%% BY COMPARING THE POSITION OF THE LIMITS WITH THE VELOCITY MINIMUM
%	%if(length(possible) > 1)
%	%	minv = find(vv(1,ilo,i) == min(vv(1,ilo,352)) );
%	%	for m = 1:length(possible)
%	%		dif(m) = possible(m) - minv;
%	%	end;
%	%end;
%	
%	minv = find(vv(1,ilo,i) == min(vv(1,ilo,i)));
%	%% TO PLOT MIN V
%	%plot(vv(1,ilo,i))
%	%hold on; plot(ilo(minv)-ilo(1)+1,vv(1,40,i),'bo')
%	if(exist('possible'));
%	if(length(possible) >= 1)
%		for m = 1:length(possible)
%			if(minv > el);
%				el = wl(end);
%				wl = wl(possible(m)+1);
%			end;
%		end;
%	end;
%	end;
%	%%%%  THERE ARE CASES WERE THERE IS ONLY ONE NEGATIVE BULGE, THIS SOLVES IT
%	if(cont == 0)
%		el = wl(end);
%	end;
%
%
%	%% EASTERN LIMIT FOUND ON TRANSPORT
%	%el = find(vtransp(ilo,i) >= 0);
%	%el = (el>wl)
%	
%	%% RECORDING ON VARIABLE
%	bczl(i,1) = wl(1);
%	bczl(i,2) = el(1);
%	%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%
%	
%	%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%	%% FIND BC VERTICAL LIMITS (depth limit)
%	dl(:,:,i) = (vv(:,:,i) > 0);
%	%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%	
%	%%%% clearing variables
%	clear possible minv;
%
%end;
%
%figure; 
%contourf(squeeze(vv(1,:,:)),10);
%hold on;
%plot(ilo(bczl(:,1)),'r-*');
%hold off;
%
%%%%% THIS LOOP WILL TAKE THE INFORMATION ABOVE, ABOUT THE
%%%%% ZONAL LIMITS OF THE BC, AND WILL CALCULATE THE TRANSPORT
%
%%% defining variable to save conputing time
%bct = 0.*ones(size(squeeze(vv(1,:,:))));
%
%%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%for i = 1:length(date)
%	for j = 1:length(lon)
%
%		%% BC depth limit
%		bcdl = find(dl(:,j,i) == 1);
%
%		if(sum(bcdl) >= 1 & size(bcdl) > 0 & lon(j) >= bczl(i,1) & lon(j) <= bczl(i,2));
%			disp('YEAH FUCKA')
%			bcdl = bcdl(end);
%			bct(j,i) = nansum(vv(1:bcdl,j,i));
%		        bct(j,i) = bct(j,i).*dep(bcdl);
%		end;
%
%	end;
%end;
%%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  IT WORKS DOWN FROM HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.
%--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.




