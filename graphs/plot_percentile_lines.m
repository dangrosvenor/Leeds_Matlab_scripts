%draw percentile lines on a frequency distribution

 prctiles=prctile(X,[10 25 50 75 90]);
 
 imed=3;
 i25=2;
 i75=4;
 
 i10=1;
 i90=5;
 
 ylims=get(gca,'ylim');
 
 iy=isnan(ydat(1).y);
 iy=find(iy==1);
 xx=xdat(1).x;
 xx(iy)='';
 
 yy=ydat(1).y;
 yy(iy)='';
 
 %median
 iprc=imed;
 h = line([prctiles(iprc) prctiles(iprc)],[ylims(1) interp1(xx,yy,prctiles(iprc))]);
 set(h,'linewidth',2);
 
 iprc=i25;
 h = line([prctiles(iprc) prctiles(iprc)],[ylims(1) interp1(xx,yy,prctiles(iprc))]);
 set(h,'linestyle','--');
 
 iprc=i75;
 h = line([prctiles(iprc) prctiles(iprc)],[ylims(1) interp1(xx,yy,prctiles(iprc))]);
 set(h,'linestyle','--');
 
 
 iprc=i10;
 h = line([prctiles(iprc) prctiles(iprc)],[ylims(1) interp1(xx,yy,prctiles(iprc))]);
 set(h,'linestyle','--');
  set(h,'color','k');
 
 iprc=i90;
 h = line([prctiles(iprc) prctiles(iprc)],[ylims(1) interp1(xx,yy,prctiles(iprc))]);
  set(h,'linestyle','--');
 set(h,'color','k');
 
 