%f=[0.001 0.01 0.1 0.5 0.99 1];

%methods={'UM pre bug-fix','UM','Robin','Robin2','Dan Robin','Dan Robin2','lognormal dN/lnD ratio','midway all or nothing'}';

col={'b','b','k','k','r','r','c','g'};
mark={'o','s','x','s','^','v','s','s'};


figure

for im=1:length(methods)
    method = methods{im};

    plot(r_ratio(end,:),r1_new(im,:),'-','color',col{im},'marker',mark{im});
    %plot(r_ratio(im,:),dm1(im,:),'-','color',col{im},'marker',mark{im});
    hold on

end

%set(gca,'xscale','log');
legend(methods);
xlabel('(rm-r1)/(r2-r1)');
ylabel('Accumulation mode size (m)');

%plot(r_ratio(end,:),r1(im,:).*ones(size(r_ratio(im,:))),'k--');
%plot(r_ratio(im,:),r2(im,:).*ones(size(r_ratio(im,:))),'k--');


figure

for im=1:length(methods)
    method = methods{im};

    plot(r_ratio(end,:),r2_new(im,:),'-','color',col{im},'marker',mark{im});
    hold on

end

plot(r_ratio(end,:),r1(end,:).*ones(size(r_ratio(end,:))),'k--');
plot(r_ratio(end,:),r2(end,:).*ones(size(r_ratio(end,:))),'k--');

set(gca,'xscale','log');
legend(methods);
xlabel('(rm-r1)/(r2-r1)');
ylabel('Coarse mode size (m)');

%%
figure

for im=1:length(methods)
    method = methods{im};

    plot(r_ratio(end,:),dn_ratio_accum(im,:),'-','color',col{im},'marker',mark{im});
    hold on

end

%plot(r_ratio(end,:),r1(end,:).*ones(size(r_ratio(end,:))),'k--');
%plot(r_ratio(end,:),r2(end,:).*ones(size(r_ratio(end,:))),'k--');

set(gca,'xscale','log');
legend(methods);
xlabel('(rm-r1)/(r2-r1)');
ylabel('dn1/(dn1+dn2)');


%%
figure

for im=1:length(methods)
    method = methods{im};

    plot(r_ratio(end,:),dm1(im,:),'-','color',col{im},'marker',mark{im});
    hold on

end

%plot(r_ratio(end,:),r1(end,:).*ones(size(r_ratio(end,:))),'k--');
%plot(r_ratio(end,:),r2(end,:).*ones(size(r_ratio(end,:))),'k--');

set(gca,'xscale','log');
legend(methods);
xlabel('(rm-r1)/(r2-r1)');
ylabel('dm1');


%%


