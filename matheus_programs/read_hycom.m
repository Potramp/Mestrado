% Read the HYCOM data processed by Edmo Campos. The program reads
% the weekly data that are written in netcdf format

% Standard depths
load stan_depth
dep=dep(1:32); % use only up to 5000 m

% HYCOM directory
prefix='/Volumes/titan/edmo/hycom/GLBa0.08/expt_18.3/plot_ferret/data_m/';

% all 3di files. Use files from 1995 to 2015.
d=dir([prefix,'archm*3di.nc']); 
id1=1531; % index for first day of 1995
idn=2753; % index for last day of 2014

% define big matrices
tem=nan*ones(32,176,idn-id1+1);
sal=tem;
%den=tem;
%uu=tem;
vv=tem;

% Define lat and lon from one netcdf file (section at -28S)
filename=[prefix,d(id1).name];
lat=ncread(filename,'Latitude');
lon=ncread(filename,'Longitude');
ilo=find(lon<=-38&lon>=-52);
ila=find(lat<=-29.98&lat>=-30.05);
start=[ilo(1) ila(1) 1 1];
count=[length(ilo) length(ila) 32 1];

for ii=id1:idn % from 1995 to 2014
filename=[prefix,d(ii).name];
disp(filename)

% Read the time from file name
yy(ii-id1+1)=str2num(filename(71:74));
gg=str2num(filename(76:78));
[a,b]=gregorian(yy(ii-id1+1),gg);
mm(ii-id1+1)=a;
dd(ii-id1+1)=b;

% Read the netcdf data 
tt=squeeze(ncread(filename,'layer_temperature',start,count))';
ss=squeeze(ncread(filename,'layer_salinity',start,count))';
%rr=squeeze(ncread(filename,'layer_density',start,count))';
%u=squeeze(ncread(filename,'u_velocity',start,count))';
v=squeeze(ncread(filename,'v_velocity',start,count))';
zz=squeeze(ncread(filename,'interface_depth',start,count))';

[m,n]=size(tt);
% Define the variables for interpolation
tti=nan*ones(32,176);
ssi=tti;
%rri=tti;
%ui=tti;
vi=tti;

for j=1:n
L=~isnan(tt(:,j));
if sum(L)~=0;
tti(:,j)=interp1(zz(L,j),tt(L,j),dep');
ssi(:,j)=interp1(zz(L,j),ss(L,j),dep');
%rri(:,j)=interp1(zz(L,j),rr(L,j),dep');
%ui(:,j)=interp1(zz(L,j),u(L,j),dep');
vi(:,j)=interp1(zz(L,j),v(L,j),dep');
end
end

% Repeat the first depth from original variables
tti(1,:)=tt(1,:);
ssi(1,:)=ss(1,:);
%rri(1,:)=rr(1,:);
%ui(1,:)=u(1,:);
vi(1,:)=v(1,:);

% Include each day on the big matrix
tem(:,:,ii-id1+1)=tti;
sal(:,:,ii-id1+1)=ssi;
%den(:,:,ii-id1+1)=rri;
%uu(:,:,ii-id1+1)=ui;
vv(:,:,ii-id1+1)=vi;
end

lat=lat(ila);
lon=lon(ilo);

save hycom_30S yy mm dd tem sal  vv lat lon dep
