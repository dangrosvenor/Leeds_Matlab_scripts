acc=0.5; %orig=1
t=f*sum(TwoDDan(idir).Q(:,:,1:6),3);

%torig=GridDan(idir).OLQBAR(:,1); %original vapour mixing ratio - Grid.OLQBAR changes with time
 tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
potemp=squeeze(sum(TwoDDan(idir).TH1,3))+tref; %potemp
            
            
i1=findheight(GridDan(idir).Z,10e3);
i2=findheight(GridDan(idir).Z,25e3);

pot=[tref(i1,1):0.25:tref(i2,1)];
            
            
%pot=[300:0.5:450];



istore=[];
for ip=2:length(pot)
%    [xx2,zz2,i]=potempsmooth(TwoDDan(idir),GridDan(idir),pot(ip),2,1);
   % i=find(potemp>pot(ip)-acc & potemp<pot(ip)+acc);
   
  if length(i)>0
   i=find(potemp>=pot(ip-1) & potemp<pot(ip));
    iorig=findheight(tref(:,1),pot(ip));
    %t_iorig=torig(idir).t(iorig);
    tot=t(i);
    %med=median(tot);
    medacc=0.05; %orig=0.1
   % med=t_iorig;
    med=prctile(t(i),20);
    d=find(abs(tot-med)/med<medacc);
    istore=[istore i(d)'];
  end
  
end
    
istore=unique(istore);

[zzi xxi]=ind2sub(size(TwoDDan(idir).TH1),istore);

xx2=GridDan(idir).Y1(xxi)/1000;
zz2=0.62+GridDan(idir).Z(zzi)/1000;
fprintf(1,'done');