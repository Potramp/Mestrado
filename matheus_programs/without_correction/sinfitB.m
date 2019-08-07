function coef=sinfitB(t,x,T)
% function coef=sinfitB(t,x,T) fits a sinusoidal curve with 
% KNOWN PERIOD to a time series. 
% x(t)=dependent variable, t=time, T=period
% newx(t)=coef(1)+coef(2)*t+coef(3)*sin((2*pi/T)*t)+coef(4)*cos((2*pi/T)*t) OR
xx=x(:);tt=t(:);
%

wt=2*pi/T;

A=[ones(size(tt)) ones(size(tt)).*tt sin(wt*tt) cos(wt*tt)];

% in sinfitb, A=[ones(size(tt)) sin(wt*tt) cos(wt*tt)];

coef=(A'*A)\(A'*xx);
