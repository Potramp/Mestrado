%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T-TIL IS MONTHLY. INTERPOLATE TO AN ANNUAL CASE
%
%
%
%
%
%
%
%
% - CORTEZIi, 2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
more('off')

%%%% LOADING T_TILDE
load ~/Dropbox/Mestrado/matheus_programs/T_tilde

rT = [T_tilde T_tilde T_tilde];
days = [1:365];
mday = [31 59 90 120 151 181 212 243 273 304 334 365];
rm = [mday mday+365 mday+365*2];
rd = [days days+365 days+365*2];

%%%% loop to make the same for all depths;
for i = 1:11
	T_tilde_y(i,:) = spline(rm,rT(i,:),rd);
end;

T_tilde_y = T_tilde_y(:,366:end-365);

save T_tilde_y.mat T_tilde_y;
