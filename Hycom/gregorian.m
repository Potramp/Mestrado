function  [m,d]=gregorian(y,jd);
% function  [m,d]=gregorian(y,jd);
% convert Julian Day to gregorian date (Month, Day) , provided the Year.
% Check also JULIAN, DAY2S and S2DAY
%
y=y(:);jd=jd(:);				% prepare the input arrays
m=zeros(1,length(y));d=m;			% preload the output arrays
% ns is the julian day of the end of each month (also for leap years)
ns=meshgrid([ 31 59 90 120 151 181 212 243 273 304 334 365]',1:length(y));
if(any(rem(y,4)==0)),
 l=find(rem(y,4)==0);
 ns(l,:)=meshgrid([31 60 91 121 152 182 213 244 274 305 335 366]',1:length(l));
end
%
m=(sum((ns< (meshgrid(jd,1:12)') )')+1)';	% get the month
ns=[zeros(1,length(y))' ns]';			% 
d=jd-(ns(m'+(1:13:length(y)*13)-1))';		% get the day
