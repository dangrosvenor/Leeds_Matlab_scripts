%plots bar graph of aerosol bins

multipop; %calculates aerosol distribution
combined_sdist33;

figure;

xdat=exp(logD(1).d);

xdat=Sc(1).s;
N(1).n(end+1)=0;

%ydat=N(1).n + N(2).n;
ydat=[Ntot 0]/1e6;
xdat=snew2*100; %convert to %

h=stairs(xdat,ydat);

for i=2:length(xdat)
    h2=line([xdat(i) xdat(i)],[0 ydat(i-1)]);
end

%set(gca,'xlim',[0 0.1e-6]);

set(gca,'xscale','log');
%set(gca,'yscale','log');

set(gca,'fontsize',20);

xlabel('Critical Supersaturation (%)');
ylabel('Number concentration (cm^{-3})');