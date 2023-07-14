function plotemmtimHf(temm,Temm,dat,col,savedir,nam,saveflag,clims,tlims,zlims,zflag)

if length(zlims)~=0;
	izmin=findheight(Temm,zlims(1));
	izmax=findheight(Temm,zlims(2));
end


if length(izmin)==0; izmin=1; end
if length(izmax)==0; izmax=length(Temm); end

if length(tlims)~=0;
	it1=findheight(temm,tlims(1));
	it2=findheight(temm,tlims(2));
end

if length(it1)==0; it1=1; end
if length(it2)==0; it2=length(temm); end


figure;
pcolor(temm(it1:it2),Temm(izmin:izmax),dat(izmin:izmax,it1:it2,col));
shading flat;
title(nam);


% if length(zlims)>0
%     set(gca,'ylim',[Temm(izmin) Temm(izmax)]);
% end
% 
% if length(tlims)>0
%     set(gca,'xlim',[temm(it1) temm(it2)]);
% end

if zflag==0 %if are using altitude
	set(gca,'ydir','reverse');
end

if exist('clims')
	if length(clims)==2;
        set(gca,'clim',clims);
	end
end

colorbar;

set(gcf,'name',nam);

if saveflag==1
	picname=[savedir nam];
	print(gcf,picname,'-dmeta');
	close(gcf);
end
