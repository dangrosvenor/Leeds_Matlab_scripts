clear fx pxx fxx diff pi pim pre

igrid=1;
dy=diff(GridDan(igrid).Y1(1:2));

itend=19;
cols={'r','b','k--'};
figure
hold on

iz=132;
iz=139;

for idir=1:3
    
	for i=1:itend
        [pxx(idir).p(:,i),fxx(idir).f(:,i)]=periodogram(TwoD_alltim(idir).W(iz,:,i),[],[],1/dy);
        f1=fft(TwoD_alltim(idir).W(iz,:,i));
        pim(i,:)=imag( f1(2:size(f1,2)/2) );
        pre(i,:)=real(f1(2:size(f1,2)/2) );
        pp(i,:)=pim(i,:).^2 - pre(i,:).^2;
        pp2(i,:)=f1.*conj(f1);
	end
	plot(2*pi*fxx(idir).f,mean(pxx(idir).p,2),cols{idir});
end

set(gca,'xlim',[5e-5 1.5e-3]);
%set(gca,'ylim',[10 6e4]);

set(gca,'yscale','log');
set(gca,'xscale','log');

