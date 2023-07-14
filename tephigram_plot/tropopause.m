function [t,h]=tropopause(t_interp,alt_interp)

hposdtdz=1200;

dtdz=diff(t_interp)./diff(alt_interp);
posdt=find(dtdz>0);
goflag=1;
ipos=1;
while goflag==1
    ihpos=posdt(ipos);
    hpos=alt_interp(ihpos);
	ihpos2=findheight(alt_interp,hpos+hposdtdz); %find point hposdtdz km above the positive lapse rate
    ipos2=ipos+ihpos2-ihpos; %index in posdt if all points hposdtdz km above also have positive dtdz
    if ipos2>length(posdt) %if have reached end of positive lapse rate points
        goflag=0;
    else
        diffipos=diff(posdt(ipos:ipos2)); %all diffs should be one if all points positive
        idiffgt1=find(diffipos>1);
        %match=posdt(ipos:ipos2)==[hpos:hpos2];
        if length(idiffgt1)==0 %if differences in posdt indices for a hposdtdz km jump are all one - i.e. dtdz positive for 1.2km
            goflag=0; %stop search - altitude we need is hpos
        else        
            ipos=ipos+idiffgt1(1); %start searching again from height where lapse rate not positive
        end
    end
end
%TTROP=t_interp(min(i));
%HTROP=alt_interp(min(i));
%TTROP=t_interp(itrop);
%HTROP=alt_interp(itrop);
t=t_interp(posdt(ipos))-273.15;
h=hpos;