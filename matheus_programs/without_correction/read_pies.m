% Read PIES data from NOAA: row format

load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/PIES_data/pies_sal.asc
%load /Users/apple/Dropbox/Mestrado/chris_programs/PIES_data/pies_sal.asc
data=pies_sal;

yy=data(2:end,1);
mm=data(2:end,2);
dd=data(2:end,3);
hh=data(2:end,4);
jd=datenum(yy,mm,dd,hh,0,0);
sal=data(2:end,5:end);
dep=data(1,5:end);

load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/PIES_data/pies_temp.asc
%load /Users/apple/Dropbox/Mestrado/chris_programs/PIES_data/pies_temp.asc
data=pies_temp;
tem=data(2:end,5:end);

clear data pies*

save piesD_ts yy mm dd hh jd sal tem dep

