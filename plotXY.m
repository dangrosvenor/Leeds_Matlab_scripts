%plots 2-D data as line plots with markers 

function [H1]=plotXY(xdat,ydat,labs,nmark) %xdat(1:n).x, ydat(1:n).y & labs(1:n).l nmark=no markers


fsize=25;


cdan(1).c=[1 0 0];  %red
cdan(3).c=[0 0.5 0.1]; %dark green
cdan(2).c=[0 0 1];  %blue
cdan(4).c=[0.5 0.5 0.5];  %grey
cdan(5).c=[0 0 0];
cdan(6).c=[0.7 0.8 0]; %yellowish

size(cdan);
scdan=ans(2);

pdan(1).p='-';
pdan(2).p=':';
%pdan(3).p='--';
%pdan(4).p='-.';

markers(1).m='h';
markers(2).m='+';
markers(3).m='o';
markers(4).m='*';
markers(5).m='.';
markers(6).m='x';
markers(7).m='s';
markers(8).m='d';
markers(9).m='^';
markers(10).m='<';
markers(11).m='>';
markers(12).m='p';
markers(13).m='h';
smark=length(markers);
spdan=length(pdan);

hold on;

ax=axes;

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
     set(H1(j).h,'LineStyle',patt(j).p);
     set(H1(j).h,'color',colour(j).c);
     

     
     markpnts=[1:round(length(xdat(j).x)/nmark):length(xdat(j).x)];
     
     H2(j).h=plot(xdat(j).x(markpnts),ydat(j).y(markpnts));
     set(H2(j).h,'LineStyle','none');
     set(H2(j).h,'color',colour(j).c);
     set(H2(j).h,'marker',mark(j).m);
     if strcmp(mark(j).m,'.')==1
         set(H2(j).h,'markersize',20);
     end
     

     
         
end

[lh ha hb hc]=legend([H1.h],labs.l,2);
for j=1:length(xdat)  
     set(ha(2*j+1),'marker',mark(j).m);
     set(ha(2*j+1),'color',colour(j).c);
     if strcmp(mark(j).m,'.')==1
        set(ha(2*j+1),'markersize',20);
     end
     
 end
 
 set(lh,'fontsize',fsize-5);
 set(gca,'fontsize',fsize);
 
