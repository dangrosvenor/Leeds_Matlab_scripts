clear diff

idir=1;

[ih ih2]=findheight(GridDan(idir).Z,0.01e3,2.5e3);
w=TwoDDan(1).W(ih:ih2,:);

for k=1:ih2-ih+1
	iwa=find(w(k,:)>1);
	wid(k)=length(iwa);
    wme(k)=mean(w(k,iwa));
end

up_wid=mean(wid) * diff(GridDan(idir).Y1(1:2)) / 1000;
mean_w=mean(wme);
    

for idir=1:length(icediagsALL)
    rho=repmat(GridDan(idir).RHON(ih:ih2),[1 size(icediagsALL(idir).i,2)]);
    
	area=icediagsALL(idir).i(ih:ih2,:,[283])/npess2(idir); %W>1_A
	aream2=diff(GridDan(1).Y1([1 end]))/1000 * area; %W>1_A
	warea_diag(idir).w=mean(area);
	
	area(find(area==0))=1;
	wmean_diag(idir).w=mean(rho.*icediagsALL(idir).i(ih:ih2,:,303)./area /npess2(idir) ); %W>1_W
    
    wflux_diag(idir).w= mean(aream2 .* rho.*icediagsALL(idir).i(ih:ih2,:,303)./area /npess2(idir) ); %W>1_W

end
