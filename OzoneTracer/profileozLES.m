scrsz=get(0,'ScreenSize');
posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

secax=1;

hf=figure('position',posit);

f=1e9*28.97/48; %converts mr into ppbv

clear cbtick
xdat(1).x=log10(max(TwoDDan(2).Q(:,:,14),[],2));
xdat(2).x=log10(max(TwoDDan(5).Q(:,:,14),[],2));

ydat(1).y=Grid.Z/1000;
ydat(2).y=Grid.Z/1000;

labs(1).l='Low updraught case';
labs(2).l='High updraught case';


is=1;
xdat(3).x=log10(data(is).dmi(9,:));
ydat(3).y=data(is).dmi(1,:)/1000;
labs(3).l='DMI flight 24th Feb';

ax=axes;
plotXY(xdat,ydat,labs,20);

cblim=get(ax,'xlim');
%cblim=f*10.^(cblim);
nn=9; %number of ticks on x-axis
xn=[cblim(1):(cblim(2)-cblim(1))/nn:cblim(2)]; %equal spacing along x-axis
xn=log10(f*10.^xn);

rd=floor(xn); %order of magnitude of values for xaxis 

for i=1:length(rd)
        xnew(i)=round(10^rd(i).*round(10.^xn(i)./10^rd(i))); %round to order of magnitude of value
end

% tvals=get(ax,'xtick');
% tvals=10.^(tvals)*f;
% 

    tvals=xnew;
    cbvals2=log10(xnew/f);

    for icb=1:length(tvals)
        te=num2str(tvals(icb),'%1.0f');
        cbtick(icb,1:length(te))=te;
    end
    
    %cbtick(1,1)=' ';
            set(ax,'xtick',cbvals2);
            set(ax,'xticklabel',cbtick);
            
            ylabel('Height (km)');
            xlabel('Max Ozone Mixing Ratio (ppbv)');

if secax==1
	% create 2nd axis for displaying dividers between ticks
	clear cb
	ax2=axes;
	%set(ax2,'ytick',[])
	set(ax2,'Color','none','xticklabel',[]); %color=none makes axis backplane transparent
	set(ax2,'yticklabel',[]); %color=none makes axis backplane transparent
	
	
	n=5; %no. of dividers between ticks
	for i=1:length(xnew)-1
        for j=1:n
            ind=(i-1)*n+j;
            cb(ind)=(xnew(i+1)-xnew(i))*j/n+xnew(i);
        end
	end
	cbnew=( (log10(cb/f)) - cblim(1) ) ./ (cblim(2)-cblim(1));
	
	set(ax2,'xtick',cbnew,'xticklabel',[],'ticklength',[0.005 0.005]); %ticklength - 1st element for 2-D, 2nd for 3-D
	
	
	xnew=get(ax,'ytick');
	
	n=5; %no. of dividers between ticks
	for i=1:length(xnew)-1
        for j=1:n
            ind=(i-1)*n+j;
            cb(ind)=(xnew(i+1)-xnew(i))*j/n+xnew(i);
        end
	end
	
	cblim=[xnew(1) xnew(end)];
	cbnew=( cb - cblim(1) ) ./ (cblim(2)-cblim(1));
	set(ax2,'ytick',cbnew,'yticklabel',[],'ticklength',[0.005 0.005]); %ticklength - 1st element for 2-D, 2nd for 3-D

end

