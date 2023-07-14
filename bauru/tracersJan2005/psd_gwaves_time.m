clear fx px diff pxx fxx

idirs=1;


itend=19;
cols={'r','b','k--'};
figure
hold on

iz=132;
for ispec=1:length(idirs)
    dy=300; %sampling period = 300 secs
	for i=1:length(GridDan(idirs(ispec)).Y1)
        [pxx(ispec).p(:,i),fxx(ispec).f(:,i)]=periodogram(TwoD_alltim(ispec).W(iz,i,1:itend),[],[],1/dy);
	end
	plot(1e3*2*pi*fxx(ispec).f,mean(pxx(ispec).p,2),cols{ispec}); %multiplied by 1e3 for easy comparison with paper
end

set(gca,'xlim',[0 15]);
%set(gca,'ylim',[10 6e4]);

set(gca,'yscale','log');
%set(gca,'xscale','log');

