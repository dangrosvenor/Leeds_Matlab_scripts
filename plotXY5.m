%plots 2-D data as line plots with markers 

function [H1,ax1,ax2]=plotXY3(xdat,ydat,labs,nmark,lwidth,leglor,logflag,xlab,ylab,ylims,zline,dual,secyA,secyB,lab2,...
    fsize,ixlab,ixdir,iydir) %xdat(1:n).x, ydat(1:n).y & labs(1:n).l nmark=no markers

%logflag=1 for log10 xaxis, =2 for yaxis & =12 for both 
%put nmark=-1 for all markers

nonew=0;


cdan(1).c=[1 0 0];  %red
cdan(3).c=[0 0.5 0.1]; %dark green
cdan(2).c=[0 0 1];  %blue
cdan(4).c=[0.5 0.5 0.5];  %grey
cdan(5).c=[0 0 0];
cdan(6).c=[0.7 0.8 0]; %yellowish
cdan(7).c=[1.0 0.7 0.7]; %salmon pink
cdan(8).c=[0.7 0.7 1.0]; %light blue
cdan(9).c=[0.6 0.6 0.8]; %turqoise?

size(cdan);
scdan=ans(2);

pdan(1).p='-';
pdan(2).p=':';
pdan(3).p='--';
pdan(4).p='-.';

% markers(1).m='d';
% markers(2).m='+';
% markers(3).m='o';
% markers(4).m='*';
% markers(5).m='.';
% markers(6).m='x';
% markers(7).m='s';
% markers(8).m='d';
% markers(9).m='^';
% markers(10).m='<';
% markers(11).m='>';
% markers(12).m='p';
% markers(13).m='h';

markers(1).m='d';
markers(2).m='o';
markers(3).m='s';


smark=length(markers);
spdan=length(pdan);

for i=1:length(labs) %remove underscores from labels as makes the next character subscript
	iund=findstr(labs(i).l,'_');
	if length(iund)>0; labs(i).l(iund)='-'; end;
end

for j=1:length(xdat)
    
            if rem(j,scdan)==0
                colour(j).c=cdan(scdan).c;
            else
                colour(j).c=cdan(rem(j,scdan)).c;
            end
            if rem(j,spdan)==0
                patt(j).p=pdan(spdan).p;
            else
                patt(j).p=pdan(rem(j,spdan)).p;
            end
            if rem(j,smark)==0
                mark(j).m=markers(smark).m;
            else
                mark(j).m=markers(rem(j,smark)).m;
            end
            
            if rem(j,3)==0
                linet(j)=3;
            else
                linet(j)=rem(j,3);
            end
            

         H1(j).h=plot(xdat(j).x,ydat(j).y);
 
        if ixdir==-1;
            set(gca,'xdir','reverse');
        end
        if iydir==-1;
            set(gca,'ydir','reverse');
        end
        
     
     hold on;
     
     if zline==1 %plot a zero line
         z=zeros([1 length(xdat(1).x)]);
         plot(z,ydat(1).y,'k');
     end
     
     set(H1(j).h,'LineStyle',patt(j).p);
     set(H1(j).h,'color',colour(j).c);
     

     if j>length(lwidth);lw=lwidth(1);else; lw=lwidth(j); end;
     set(H1(j).h,'linewidth',lw);     

     if nmark==-1
         markpnts=[1:length(xdat(j).x)];
     elseif nmark==0
         markpnts=[];
     else
        markpnts=[1:round(length(xdat(j).x)/nmark):length(xdat(j).x)];
     end
    
         H2(j).h=plot(xdat(j).x(markpnts),ydat(j).y(markpnts));

     
     set(H2(j).h,'LineStyle','none');
     set(H2(j).h,'color',colour(j).c);
     set(H2(j).h,'marker',mark(j).m);
%     if strcmp(mark(j).m,'.')==1
         set(H2(j).h,'markersize',5);
         set(H2(j).h,'markerfacecolor',colour(j).c);
 %    end
     

      
         
end

if length(ylims)>1
     set(gca,'ylim',ylims);
end
 

if dual==1
 
 ynew=get(gca,'ytick');
 xnew=get(gca,'xtick');
 ax1=gca;
 
 if logflag==2 | logflag==12
     set(ax1,'yscale','log');
 else
    if nonew==0;
        posgca=get(gca,'position');
        ax2=axes('position',posgca);
 %       set(ax2,'Color','none','xticklabel',[],'xtick',0,'box','off'); %color=none makes axis backplane transparent
        set(ax2,'xticklabel',[],'xtick',0,'box','off'); %color=none makes axis backplane transparent

    end
    
    set(ax2,'yticklabel',[]); %color=none makes axis backplane transparent

    clear cbnew cb

    cb=ynew;
	
	cblim=get(ax1,'ylim');
	cbnew=( cb - cblim(1) ) ./ (cblim(2)-cblim(1));
    
    for i=1:length(ynew)
        ip=findheight(secyA,ynew(i));
        if ynew(i)<=secyA(end)
            te=num2str(round2(secyB(ip),0));
            ctickstr(i,1:length(te))=te;
        else
            ctickstr(i,1)=' ';
        end
    end
    
    
	set(ax2,'ytick',cbnew,'yticklabel',ctickstr); %ticklength - 1st element for 2-D, 2nd for 3-D
    set(ax2,'yaxislocation','r','yminortick','off');
    
    set(ax2,'fontsize',fsize);
    ylabel(lab2);
  end
  
else
    if logflag==2 | logflag==12
         set(ax1,'yscale','log');
    end
end



if logflag==1 | logflag==12
  set(gca,'xscale','log');
end

if logflag==2 | logflag==12
  set(gca,'yscale','log');
end


if exist('ax2')
    %switch_objects_depth(gcf,ax1,lh); %switchs order of figure children (ax and ax2) so that can zoom on ax
    axes(ax1);
    set(ax1,'Color','none'); %makes backpane of axis one transparent
 
    %axes(lh);
    %switch_objects_depth(gcf,ax2,lh);
    %set(ax1,'color','none'); %color=none makes axis backplane transparent
    %


end

set(gca,'fontsize',fsize);

[lh ha hb hc]=legend([H1.h],labs.l,leglor);
for j=1:length(xdat)
     if nmark~=0
        set(ha(2*j+1),'marker',mark(j).m);
     end
     set(ha(2*j+1),'color',colour(j).c);
     if strcmp(mark(j).m,'.')==1
        set(ha(2*j+1),'markersize',20);
     end
     
 end
 
 set(lh,'fontsize',fsize+2);
 
 if ixlab==1
     xlabel(xlab);
 end
 ylabel(ylab);
 
 

 
%  if logflag==1 | logflag==12
%      set(gca,'xscale','log');
%  else %extra ticks on x-axis
%      if nonew==0
%     ax2=axes;
%     set(ax2,'xticklabel',[],'ytick',0); 
%     set(ax2,'yticklabel',[]); %color=none makes axis backplane transparent
% 
%     clear cbnew cb
% 	
% 	n=5; %no. of dividers between ticks
% 	for i=1:length(xnew)-1
%         for j=1:n
%             ind=(i-1)*n+j;
%             cb(ind)=(xnew(i+1)-xnew(i))*j/n+xnew(i);
%         end
% 	end
% 	
% 	cblim=[xnew(1) xnew(end)];
% 	cbnew=( cb - cblim(1) ) ./ (cblim(2)-cblim(1));
% 	set(ax2,'xtick',cbnew,'xticklabel',[],'ticklength',[0.005 0.005]); %ticklength - 1st element for 2-D, 2nd for 3-D
%     
%     nonew=1;
% end
% end


  
