%plots 2-D data as line plots with markers

function [H1,ax1,ax2]=plotXY6(xdat,ydat,labs,nmark_array,lwidth,leglor,logflag,xlab,ylab,ylims,zline,dual,secyA,secyB,lab2,...
    fsize,ixlab,ixdir,iydir,xloc,time_highlight_path,highlight_type,ierror_bars,errordatU,errordatL,marksize,ichoose_styles,pdan,cdan,mdan,line_widths,iovr_leg_line) %xdat(1:n).x, ydat(1:n).y & labs(1:n).l nmark=no markers

%logflag=1 for log10 xaxis, =2 for yaxis & =12 for both
%put nmark=-1 for all markers
%has some provision for dual axes (dual=2). Think need to set xloc to say
%which axis want to plot which xdat on. E.g. xloc=[1 1 0 0]; But not tested
%recently (as of July, 2016)


error_bars=0;

%calculate the total span of xdat for plotting error bars of a constant
%width
minX=9e99; maxX=-9e99;
for idat=1:length(xdat)
    minX = min([min(xdat(idat).x) minX]);
    maxX = max([max(xdat(idat).x) minX]);
end
spanX = maxX-minX;


if length(time_highlight_path)>1
    
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
    h=(ylims(2)-y)-dy*sc_f;


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

if ichoose_styles==0

    cdan(1).c=[1 0 0];  %red
    cdan(3).c=[0 0.5 0.1]; %dark green
    cdan(2).c=[0 0 1];  %blue
    cdan(4).c=[0.5 0.5 0.5];  %grey
    cdan(5).c=[0 0 0];
    cdan(6).c=[0.7 0.8 0]; %yellowish
    cdan(7).c=[1.0 0.7 0.7]; %salmon pink
    cdan(8).c=[0.7 0.7 1.0]; %light blue
    cdan(9).c=[0.6 0.6 0.8]; %turqoise?



    if lwidth==0
        pdan(1).p='none';
        pdan(2).p='none';
        pdan(3).p='none';
        lwidth=1;
    else
        pdan(1).p='-';
        pdan(2).p='--';
        pdan(3).p='-.';
        %pdan(4).p='-.';
    end





    
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

size(cdan);
scdan=ans(2);

smark=length(markers);
spdan=length(pdan);



end

if lwidth==0
    pdan(1).p='none';
    pdan(2).p='none';
    pdan(3).p='none';
    lwidth=1;    
end



leg_length_tot=0;
for i=1:length(labs) %remove underscores from labels as makes the next character subscript
	iund=findstr(labs(i).l,'_');
%	if length(iund)>0; labs(i).l(iund)='-'; end;
    leg_length_tot = leg_length_tot + length(labs(i).l);
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
%        set(ax2,'xaxislocation','t');
        set(ax2,'yaxislocation','r');
        
        set(ax2,'fontsize',fsize);
%        ylabel(lab2);

          
 
end

        
for j=1:length(xdat)
    
        nmark=nmark_array(j);
    
        if ichoose_styles==0
            
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
            
        else
            colour(j).c=cdan(j).c;
            patt(j).p=pdan(j).p;
            mark(j).m=mdan(j).m;
            if isnan(line_widths(1).l)==0 %set at NaN at start of watervap
                lwidth(j)=line_widths(j).l;  %so, only change if was overwritten in watervap
            end
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
     
% **** Plot commands and error bars   *****************************     


         H1(j).h=plot(xdat(j).x,ydat(j).y);
         

         
         switch ierror_bars
             case 'horiz'

                 hold on

                 Herr(j).h=herrorbar(xdat(j).x,ydat(j).y,errordatU(j).dat,errordatL(j).dat);

                 % adjust horizontal error bars to remove the "tails" at the end
                 errorbarYData          = ...
                     get(Herr(j).h, 'YData'       );
                 
                 errorbarXData          = ...
                     get(Herr(j).h, 'XData'       );

                 errorbarYData = errorbarYData{1};
                 

                 errorbarYData(1:9:end) = ...
                     errorbarYData(4:9:end);
                 errorbarYData(2:9:end) = ....
                     errorbarYData(4:9:end);
                 errorbarYData(7:9:end) = ...
                     errorbarYData(4:9:end);
                 errorbarYData(8:9:end) = ...
                     errorbarYData(4:9:end);
                 set(Herr(j).h, 'YData', errorbarYData,'XData',errorbarXData{1});

                 error_bars=1;
                 
             case 'horiz2'
                 %Daniel McCoy;s script
                 Herr(j).h = errorbarYY('horiz',xdat(j).x,ydat(j).y,errordatU(j).dat,gca,colour(j).c,'none',1,0.005);
                 error_bars=1;

             case 'vert'
                 hold on

                 Herr(j).h=errorbar(xdat(j).x,ydat(j).y,errordatU(j).dat,errordatL(j).dat);

                 hE_c  =  get(Herr(j).h     , 'Children'    );

                 % adjust horizontal error bars to remove the "tails" at the end
                 errorbarXData          = ...
                     get(hE_c(2), 'XData'       );

%                 errorbarXData = errorbarXData{1};


                 errorbarXData(1:9:end) = ...
                     errorbarXData(4:9:end);
                 errorbarXData(2:9:end) = ....
                     errorbarXData(4:9:end);
                 errorbarXData(7:9:end) = ...
                     errorbarXData(4:9:end);
                 errorbarXData(8:9:end) = ...
                     errorbarXData(4:9:end);
                 set(Herr(j).h, 'XData', errorbarXData);

                 error_bars=1;
                 
             case 'vert2'
                 %Daniel McCoy;s script
                 Herr(j).h = errorbarYY('vert',xdat(j).x,ydat(j).y,errordatU(j).dat,gca,colour(j).c,'none',1,0.005,spanX);
                 error_bars=1;
                 
             case 'vert2 filled'
                 %Daniel McCoy;s script
                 hold on
                 Herr(j).h = errorbarYY('vert filled bar',xdat(j).x,ydat(j).y,errordatU(j).dat,gca,colour(j).c,'none',1,0.005,spanX);
                 error_bars=0;                 

         end
         


% -----------------------------------------------------

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
     
     if error_bars==1
%         set(Herr(j).h,'LineStyle',patt(j).p);
%         switch ierror_bars
%             case 'horiz2'
%                 for im=1:length(Herr(j).h)
%                     set(Herr(j).h(im),'color',colour(j).c);
%                     set(Herr(j).h(im),'linewidth',lw/2); 
%                 end
%             otherwise
                 set(Herr(j).h,'color',colour(j).c);
                 set(Herr(j).h,'linewidth',lw/2); 
%         end
     end


  
    
     

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
         set(H2(j).h,'markersize',marksize);  %for the sctual markers
         set(H2(j).h,'markerfacecolor',colour(j).c);
         set(H2(j).h,'markeredgecolor','k');         
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
            ip=findheight_nearest(secyA,ynew(i));
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
        set(ax2,'yminortick','off');
        set(ax1,'yminortick','off');
        
    end
 
    %axes(lh);
    %switch_objects_depth(gcf,ax2,lh);
    %set(ax1,'color','none'); %color=none makes axis backplane transparent
    %


end

set(gca,'fontsize',fsize);

%now for the legend
if leg_length_tot>0 & leglor~=-99 %was there any writing in the legend? - if not then don't plot. Or if leg swithed off
    ileg_count=0;
    for ileg_loop=1:length(H1)  %
       if isnan(labs(ileg_loop).l)==0
           ileg_count=ileg_count+1;
           labs_leg(ileg_count).l = labs(ileg_loop).l;
           H1_leg(ileg_count).h = H1(ileg_loop).h;       
       end
    end
%    [lh ha hb hc]=legend([H1.h],labs.l,leglor); %write legend

    [lh ha hb hc]=legend([H1_leg.h],labs_leg.l,leglor); %write legend    
    set(lh,'fontsize',fsize-2); %need to put this before the marker labelling for the legend
end

for j=1:length(xdat)    
     nmark=nmark_array(j);
     
     %index for the markers for the legend - the middle one
         %order is length(xdat) where markers cannot be set, then
         %outer 1, mid 1, outer 2, mid 2, etc.
         ha_ind = length(xdat)+2*(j-1)+2;

         if iovr_leg_line==1
          set(ha(ha_ind-1),'linestyle','-');
         end

         
     if nmark~=0 & leg_length_tot>0 & leglor~=-99 %if there was any writing in the legend - if not then don't plot
        
         
       %  set(ha(2*j+3),'marker',mark(j).m);
                 %%set(ha(2*j+3),'marker',mark(j).m);  %put the markers on the legend line
                set(ha(ha_ind),'marker',mark(j).m);  %put the markers on the legend line                 
         %(ha(length(xdat)+2*j-1 = marker on both sides of legend line, ha(length(xdat)+2*j =
         % marker in centre of line

%         set(ha(2*j+3),'color',colour(j).c);
         set(ha(ha_ind),'color',colour(j).c);  
         set(ha(ha_ind),'markerfacecolor',colour(j).c);
%         set(ha(ha_ind),'markeredgecolor','k');  

        
         if strcmp(mark(j).m,'.')==1
             set(ha(2*j+3),'markersize',20); %for the legend
             %        set(ha(1+2*j),'markersize',20);
         end

     end
     
 end
 

 
 if ixlab==1
     xlabel(xlab,'verticalalignment','top');
 end
 ylabel(ylab);
 
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







  
