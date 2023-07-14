acc=1;
t=f*sum(TwoDDan(idir).Q(:,:,1:6),3);
pot=[300:0.5:450];

tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
potemp=squeeze(sum(TwoDDan(idir).TH1,3))+tref; %potemp
%potemp=squeeze(sum(TwoDDan(idir).Q(:,:,1),3))+tref; %potemp


istore=[];
for ip=1:length(pot)
%    [xx2,zz2,i]=potempsmooth(TwoDDan(idir),GridDan(idir),pot(ip),2,1);
    i=find(potemp>pot(ip)-acc & potemp<pot(ip)+acc);
    tot=t(i);
    med=median(tot);
    d=find(abs(tot-med)<0.1);
    istore=[istore i(d)'];
end
    
istore=unique(istore);

[zzi xxi]=ind2sub(size(TwoD.TH1),istore);

xx=GridDan(idir).Y1(xxi)/1000;
zz=0.62+GridDan(idir).Z(zzi)/1000;
fprintf(1,'done');