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

%%%% LOADING S_TILDE
load ~/Dropbox/Mestrado/matheus_programs/S_tilde

rT = [S_tilde S_tilde S_tilde];
days = [1:365];
mday = [31 59 90 120 151 181 212 243 273 304 334 365];
rm = [mday mday+365 mday+365*2];
rd = [days days+365 days+365*2];

%%%% loop to make the same for all depths;
for i = 1:11
	S_tilde_y(i,:) = spline(rm,rT(i,:),rd);
end;

S_tilde_y = S_tilde_y(:,366:end-365);

save S_tilde_y.mat S_tilde_y;
