% read woa09 data downloaded from 
%http://www.nodc.noaa.gov/
% Read the readme file.
%
% these data are distributed as ascii files, in a 180x360
% maps for 33 levels (5500m) at standard depths for T. 
% Obs.: anual data have 33 levels, the rest has 24.
%

% Get the  file names from the directory
   d=dir('*.dat');

for i=1:length(d);
str=['load ',d(i).name];
eval(str)

disp(d(i).name);

% define the files with a common name:

str=['data=',d(i).name(1:end-4),';'];
eval(str)

clear *an1

% Define which kind of variable is read: temp or salt. 
% This will be used to save the files. 

  if d(i).name(1)=='t' 
  var=['tem'];
  else
    var=['sal'];
end

% Clear the flagged data
L=data<-90;
data(L)=nan;

% To put in map format

[m,n]=size(data);
data=reshape(data',360,180,m*n/360/180);

% Change the name of the variable
data=permute(data,[2 1 3]);
data=[data(:,181:360,:) data(:,1:180,:)];
str=[var,'=data;'];
eval(str)

%lat=[-89.5:89.5];
%lon=[.5:360];

str=['save ',d(i).name(1:end-4),' ', var,';'];
eval(str)
end
