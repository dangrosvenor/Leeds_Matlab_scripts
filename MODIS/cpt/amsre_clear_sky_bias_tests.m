%Run plot_global and save the fields (P_save) for different CF thresholds

% lwp_3pct, lwp_5pct and lwp_7pct

figure
plot(lwp_7pct(:),lwp_5pct(:),'kx'); hold on

xlabel('LWP for 7% CF threshold plot (g m^{-2})');
ylabel('LWP for 5% CF threshold plot (g m^{-2})');
%One-to-one line
x=[-20:0.1:50];
y=x;
plot(x,y,'b-');

xx=lwp_7pct(:);
yy=lwp_5pct(:);
i=find(isnan(xx)==1 | isnan(yy)==1);
xx(i)=[];
yy(i)=[];
p=corr(xx(:),yy(:));
title(['Correlation coeff = ' num2str(p)]);



figure
plot(lwp_7pct(:),lwp_3pct(:),'kx'); hold on

xlabel('LWP for 7% CF threshold plot (g m^{-2})');
ylabel('LWP for 3% CF threshold plot (g m^{-2})');
%One-to-one line
x=[-20:0.1:50];
y=x;
plot(x,y,'b-');

xx=lwp_7pct(:);
yy=lwp_3pct(:);
i=find(isnan(xx)==1 | isnan(yy)==1);
xx(i)=[];
yy(i)=[];
p=corr(xx(:),yy(:));
title(['Correlation coeff = ' num2str(p)]);