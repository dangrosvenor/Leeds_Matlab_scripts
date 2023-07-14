%plots 2-D data as line plots with markers

function [H1,ax1,ax2]=plotXY6_3D(xdat,ydat,zdat,labs,nmark_array,lwidth,leglor,logflag,xlab,ylab,zlab,ylims,zline,dual,secyA,secyB,lab2,...
    fsize,ixlab,ixdir,iydir,xloc,time_highlight_path,highlight_type) %xdat(1:n).x, ydat(1:n).y & labs(1:n).l nmark=no markers

%logflag=1 for log10 xaxis, =2 for yaxis & =12 for both
%put nmark=-1 for all markers

marker_size=10*ones([1 length(xdat)]);
%marker_size(1)=15;
marker_size(end)=15;



if length(time_highlight_path)>1
    
    marker_size(end-1)=15;
    marker_size(end)=10;
    
    x=time_highlight_path(1);
    
    if length(ylims)<1
        
        ymin=9e99;
        ymax=-9e99;
        for i=1:length(ydat)
            ymin=min([ymin min(ydat(i).y)]);
            ymax=max([ymax max(ydat(i).y)]);
        end
        
        dy=ymax-ymin;
        ylims=[ymin-dy*0.1 ymax+dy*0.1];
        
        
        
    end
    
%    y=ylims(1);
%    h=ylims(2)-ylims(1);

%    y=ymin;
%    h=ymax-ymin;
    
    sc_f=0.005;
    dy=ylims(2)-ylims(1);
    y=ylims(1)+dy*sc_f;
        
%        y=ymin;
%        h=ymax-ymin;        
    h=(ylims(2)-ylims(1))-dy*sc_f-y;


    w=time_highlight_path(2)-time_highlight_path(1);


    X=[x x+w x+w x];
    Y=[y y y+h y+h];





    switch highlight_type
        case 'fill'
%                fill(X/24,Y,'c','linestyle','none','facealpha',0.3);
                H_fill=fill(X/24,Y,'c','linestyle','none');                
        case 'rectangle'
            %        rectangle('Position',[x/24,y,w/24,h],'LineWidth',2,'LineStyle','--');
            line(X([1 4])/24,Y([1 4]),'LineWidth',2,'LineStyle','--','color','k');
            line(X([2 3])/24,Y([2 3]),'LineWidth',2,'LineStyle','--','color','k');
    end
    hold on

end


nonew=0;


if length(nmark_array)<length(xdat)
  nmark_array=ones([1 length(xdat)])*nmark_array(1);
end

cdan(1).c=[1 0 0];  %red
cdan(3).c=[0 0.5 0.1]; %dark green
cdan(2).c=[0 0 1];  %blue
cdan(4).c=[0.5 0.5 0.5];  %grey
cdan(5).c=[0 0 0];
cdan(6).c=[0.7 0.8 0]; %yellowish
cdan(7).c=[1.0 0.7 0.7]; %salmon pink
cdan(8).c=[0.7 0.7 1.0]; %light blue
cdan(9).c=[0.6 0.6 0.8]; %turqoise?

%cdan(length(xdat)).c='none';

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

for idat=1:length(xdat)
	markers(idat).m='o';
end




smark=length(markers);
spdan=length(pdan);

for i=1:length(labs) %remove underscores from labels as makes the next character subscript
	iund=findstr(labs(i).l,'_');
%	if length(iund)>0; labs(i).l(iund)='-'; end;
end

if dual~=2
    dats=[1:length(xdat)];
end


if dual==2
    
 ax1=axes;
 ax2=axes;
 
% ynew=get(gca,'ytick');
% xnew=get(gca,'xtick');
% ax1=gca;

 %       set(ax2,'Color','none','xticklabel',[],'xtick',0,'box','off'); %color=none makes axis backplane transparent
 %      set(ax2,'xticklabel',[],'xtick',0,'box','off'); %color=none makes axis backplane transparent

    
 %   set(ax2,'yticklabel',[]); 

	
        
%		set(ax2,'ytick',cbnew,'yticklabel',ctickstr); %ticklength - 1st element for 2-D, 2nd for 3-D
        set(ax2,'xaxislocation','t');
        set(ax2,'yaxislocation','r');
        
        set(ax2,'fontsize',fsize);
%        ylabel(lab2);

          
 
end

        
for j=1:length(xdat)
    
        nmark=nmark_array(j);
    
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
      
     if dual==2
         if xloc(j)==1
             axes(ax1);
         else
             axes(ax2);
         end
     end

         H1(j).h=plot3(xdat(j).x,ydat(j).y,zdat(j).y);
 
        if ixdir==-1;
            set(gca,'xdir','reverse');
        end
        if iydir==-1;
            set(gca,'zdir','reverse');
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
    
         H2(j).h=plot3(xdat(j).x(markpnts),ydat(j).y(markpnts),zdat(j).y(markpnts));

     
     set(H2(j).h,'LineStyle','none');
     set(H2(j).h,'color',colour(j).c);
     set(H2(j).h,'marker',mark(j).m);
%     if strcmp(mark(j).m,'.')==1
         set(H2(j).h,'markersize',marker_size(j));
     
         
         if j~=length(xdat)+1
             set(H2(j).h,'markerfacecolor',colour(j).c);
         end
         
     
 %    end
     

      
         
end

if length(ylims)>1
    if dual==2
        axes(ax2);
        set(gca,'ylim',ylims);
        axes(ax1);
    end
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
    
    set(ax2,'yticklabel',[]); 

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
    
    if dual==2
        set(ax2,'xaxislocation','t','xminortick','off');
        
    end
 
    %axes(lh);
    %switch_objects_depth(gcf,ax2,lh);
    %set(ax1,'color','none'); %color=none makes axis backplane transparent
    %


end

set(gca,'fontsize',fsize);

[lh ha hb hc]=legend([H1.h],labs.l,leglor); %write legend
 set(lh,'fontsize',fsize+2); %need to put this before the marker labelling for the legend
for j=1:length(xdat)
     nmark=nmark_array(j);
     if nmark~=0
%        set(ha(2*j+1),'marker',mark(j).m);
        set(ha(2*j+length(xdat)),'marker',mark(j).m);  %put the markers on the legend line
        %(ha(length(xdat)+2*j-1 = marker on both sides of legend line, ha(length(xdat)+2*j = 
        % marker in centre of line

    %     set(ha(2*j+1),'color',colour(j).c);
         set(ha(2*j+length(xdat)),'color',colour(j).c);
         if strcmp(mark(j).m,'.')==1
    %        set(ha(2*j+1),'markersize',20);
            set(ha(2*j+length(xdat)),'markersize',20);
         end
     
     end
     
 end
 

 
 if ixlab==1
     xlabel(xlab);
 end
 ylabel(ylab);
 zlabel(zlab);
 
%set these if they don't exist as some matlab versions give an error if not returned
if ~exist('ax1')
	ax1=0;
end

if ~exist('ax2')
	ax2=0;
end

 
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







  
