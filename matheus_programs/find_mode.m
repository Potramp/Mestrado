function [mask2,mask1]=find_mode(var1,minval,maxval,var2);
% This is a function to mark the profiles for the definition of the mode water.
% var: variable
% lon, lat, pre: spacial axis
% minval and maxval: range for the search

[nl,nn]=size(var1);

% map of number of observations
mask1=nan*ones(size(var1));
mask2=nan*ones(size(var2));

for i=1:nl
aux1=var1(i,:);
aux2=var2(i,:);
L=aux1>=minval&aux1<=maxval;
if sum(L)>0
mask1(i,L)=aux1(L);
mask2(i,L)=aux2(L);
end
end

