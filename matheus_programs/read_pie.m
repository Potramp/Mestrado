function [yy mm dd jd dep tem sal] = read_pie(x)
% Read PIES data from NOAA: row format

letra=['ABCD'];
str=(['load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/PIES_data/pies_tem_' letra(x) '.asc']);
%str=(['load /Users/apple/Dropbox/Mestrado/chris_programs/PIES_data/pies_tem_' letra(x) '.asc']);
eval(str);
str=(['data=pies_tem_' letra(x) ';']);
eval(str);

tem=data(2:end,5:end);
yy=data(2:end,1);
mm=data(2:end,2);
dd=data(2:end,3);
hh=data(2:end,4);
jd=datenum(yy,mm,dd,hh,0,0);
sal=data(2:end,5:end);
dep=data(1,5:end);

clear data pies*

%save piesD_ts yy mm dd hh jd sal tem dep

