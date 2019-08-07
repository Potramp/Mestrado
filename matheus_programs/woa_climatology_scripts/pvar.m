function pv=pvar(x,y);
% pv is the percent of variance of x explained by y
% x is the time series with the whole signal that may include y 
% y is the time series taken as a component of x
pv=100*(1-var(x(:)-y(:))/var(x(:))); 
