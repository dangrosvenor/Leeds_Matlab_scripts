%Tperts_for_paper   *** 16th Aug, 2007 ***

iz=findheight(GridDan(1).Z,1000);

[Tpert_max itmax]=maxALL(tpert_max(1).t(1:iz,1:6)); %dump 6 is the end of heating period

iz2=findheight(GridDan(1).Z,400);
[Vpert_max ivmax]=maxALL(vappert_max(1).dat(1:iz2,1:6)); %dump 6 is the end of heating period

iz1=3;
xinds=6:9;
xinds=5:11;
%xinds=8;
f=1e6*28.97/18;

for idir=1:3

for it=1:6
    p_dat=GridDan(1).PREFN(iz1:end);
    
    if length(size(tpertTimH_full(idir).t))==4
        t_dat=mean(mean(tpertTimH_full(idir).t(iz1:end,xinds,xinds,it),2),3) + tpertTimH_full(idir).med(iz1:end,it);    
        q_dat=mean(mean(vappertTimH_full(idir).t(iz1:end,xinds,xinds,it),2),3) + vappertTimH_full(idir).med(iz1:end,it);
    else
        t_dat=mean(mean(tpertTimH_full(idir).t(iz1:end,xinds,it),2),3) + tpertTimH_full(idir).med(iz1:end,it);    
        q_dat=mean(mean(vappertTimH_full(idir).t(iz1:end,xinds,it),2),3) + vappertTimH_full(idir).med(iz1:end,it);
	end

    qsat_dat=satvappress(t_dat,'goff','liq',p_dat,1)/f;
    alt_dat=GridDan(1).Z(iz1:end);
    
	[CAPE,CIN,HLCL,TLCL,PLCL]=calc_cape(p_dat,t_dat,q_dat,qsat_dat,alt_dat); 
    
  %  plot_tephi_data2(t_dat,p_dat,q_dat,qsat_dat,alt_dat);
    
    capes(idir).cape(it)=CAPE;
    
    
end

end