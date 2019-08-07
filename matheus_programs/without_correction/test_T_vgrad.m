%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODE WATERS MUST HAVE A VERTICAL TEMPERATURE GRADIENT LESS INTENSE
% THAN THE SEASONAL AND PERMANENT THERMOCLINES. WE CALCULATE THE GRADIENT
% AND AVERAGE IT FOR EVERY PIES STATION TO DETERMINE A GOOD GRADIENT 
% PARAMETER FOR THE THERMOCLINE.
%
% -CORTEZI, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
more('off');

%%%% LOADING VARIABLES
load PV_in_vars;
load pies_latlon;
load CTD_and_Argo_p_levels;

%%%% CUT VARS BELOW 600 dbar
Pdenseries = Pdenseries(1:36,:,:);
Ts = Ts(1:36,:,:);
Ss = Ss(1:36,:,:);
dep = dep(1:36);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('1. Vertical gradient of T for thermocline...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% CALCULATING POTENTIAL TEMP.
for i = 1:4
        potTs(:,:,i) = sw_ptmp(Ss(:,:,i),Ts(:,:,i),dep',0);
end;

%%%% making a matrix the same size as the temporal series but for depths(pressure)
zz = repmat(dep',[1 2429 4]);

%%%% making depth out of pressure
zz = sw_dpth(zz,pies_lat_lon(1,1));

%%%% zz is depth, so should be negative.
T_vgrad = diff(potTs,1,1)./diff(-zz);

    
%-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('2. Average vertical gradient of T for each station...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_vgrad = squeeze(nanmean(T_vgrad,2));

%%%% DEFINING LETTERS FOR EACH STATION
letra = ['ABCD'];


for i = 1:4;
	subplot(2,2,i)
	plot(T_vgrad(:,i),-(dep(1:end-1)+dep(2:end)./2),'linewidth',2);
	if(i >= 3);
		xlabel(['Temperature vertical gradient.'])
	end;
	if(rem(i,2) ~= 0);
		ylabel('Depth');
	end;
	title(['PIES ' letra(i)]);
	set(gca,'xtick',[0:0.005:0.04])
	grid on;
	xlim([0 0.04]);ylim([-1000 0]);
end;

axis 'tight';
print -depsc2 T_vgrad_lim.eps
print -dpng T_vgrad_lim.png











