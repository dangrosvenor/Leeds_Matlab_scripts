tt=[min(t(i)):0.05:max(t(i))*1.02];
for ii=2:length(tt)
    binpts=find(t(i)>=tt(ii-1) & t(i)<tt(ii));
    n(ii)=length(binpts);
end

figure
plot(tt,n);