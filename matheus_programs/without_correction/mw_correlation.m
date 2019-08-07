%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE OUTPUT OF mw_interp IS LOADED AND CORRELATION WITH
% THE OUTPUT OF transport.m (on ~/Dropbox/Mestrado/Hycom)
% IS CALCULATED.
%
%
%
% -CORTEZI, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
more('off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. loading data...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%% LOAD OUTPUTS
%% INTERPOLATED
load ~/Dropbox/Mestrado/matheus_programs/without_correction/interp_mw_thick
%% RAW
load ~/Dropbox/Mestrado/matheus_programs/without_correction/mw_thickness_vars
%% TIMEAXIS FOR PLOTS
load ~/Dropbox/Mestrado/matheus_programs/without_correction/timeaxis
%% HYCOM TRANSPORT
load ~/Dropbox/Mestrado/Hycom/bc_transp    
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Correlation...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:4
	junk(:,i) = corr(thick(:,i), V);
end;

