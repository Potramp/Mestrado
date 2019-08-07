function [yy mm dd jd dep sal ] = S_read_pie(x)
% Read PIES data from NOAA: row format

letra=['ABCD'];
str=(['load /Users/matheusvcortezi/Dropbox/Mestrado/chris_programs/PIES_data/pies_sal_' letra(x) '.asc']);
eval(str);
str=(['data=pies_sal_' letra(x) ';']);
eval(str);

sal=data(2:end,5:end);
yy=data(2:end,1);
mm=data(2:end,2);
dd=data(2:end,3);
hh=data(2:end,4);
jd=datenum(yy,mm,dd,hh,0,0);
dep=data(1,5:end);

clear data pies*

%save piesD_ts yy mm dd hh jd sal tem dep

