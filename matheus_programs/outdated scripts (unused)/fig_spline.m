more 'off'

%levels = 100:100:3000


for i = 1:10:100;
	indx=find(pres==prange(i));
	figure; hold on;
	plot(tau(indx),temp(indx), 'ok');
	plot(taurange,tempergrid(i,:));
	xlabel('tau');
	ylabel('temp');
	title(['T and T interp for p=', num2str(prange(i))]);
	str=(['print -dpng splinep', num2str(prange(i))]);
	eval(str);

end;


