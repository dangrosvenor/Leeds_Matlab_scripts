function [minlnb,maxlnb,meanlnb_abv,meanlnb_bel,lnbbins_pos,lnbbins_neg,bins]=...
    lnb_calcs(lnb2d,GridDan,add_ground_height)

minlnb=min(lnb2d,[],2);
maxlnb=max(lnb2d,[],2);

zref=repmat(GridDan(1).Z(1:size(lnb2d,1))/1000+add_ground_height,[1 size(lnb2d,2)]);
lnbdiff=lnb2d-zref;
meanlnb_abv=zref(:,1)+meanselect(lnbdiff,'dat>0'); %calculate the mean only for points where lnb is lower than where air at
meanlnb_bel=zref(:,1)+meanselect(lnbdiff,'dat<0'); %for points where lnb is higher 

bins=GridDan(1).Z/1000+add_ground_height;
ipos=find(lnbdiff>=0);
ineg=find(lnbdiff<0);

lnbtemp=lnb2d;
lnbtemp(ineg)=0;
lnbbins_pos=binner(lnbtemp,bins)'; %put lnbs into bins - positiviely buoyant only

lnbtemp=lnb2d;
lnbtemp(ipos)=0;
lnbbins_neg=binner(lnbtemp,bins)'; %put lnbs into bins - negatively buoyant only