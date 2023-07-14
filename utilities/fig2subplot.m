function [hsub_new,Hcb2,pos_sub,pos_cb,Hleg2]=fig2subplot(hfig,hsub,fontsize,colbar_loc)

hfig_new = get(hsub,'parent');


%get the axis handle from the figure to be copied
Haxs = findobj(hfig,'type','axes');
 for i=1:length(Haxs)
        if ~isa(handle(Haxs(i)),'scribe.colorbar');
            hax=Haxs(i);
%            if isequal(double(H.axes),Ha1)
%                Hc1s=[Hc1s,Haxs(i)];
%            end
        end
 end
    
%hf2=figure(2);  %will get the user to provide the subplot handle
%s1=subplot(321);
pos_sub=get(hsub,'Position');
delete(hsub);
hsub_new=copyobj(hax,hfig_new);
set(hsub_new, 'Position', pos_sub);

%now copy the colorbar
Hcb1 = find_peer_colorbars_of_an_axes(hax);
if length(colbar_loc)>0
    set(Hcb1,'location',colbar_loc)
end
Hcb2 = copyobj(Hcb1,hfig_new);
%change its size and position to underneath the required subplot
%axis position is [left bottom width height]
pos_cb = [pos_sub(1)+pos_sub(3)*0.05 pos_sub(2)-pos_sub(4)*0.5 pos_sub(3)*0.9 pos_sub(4)*0.1];
%set(Hcb2,'position',pos_cb);

%%

%now copy the legend
Hall=findobj(hfig);
Hleg=NaN;
imark=1;
for ileg=1:length(Hall)
   type = get(Hall(ileg),'type');
   if strcmp(get(Hall(ileg),'tag'),'legend')==1
      Hleg = Hall(ileg);
   end
   if strcmp(type,'line')==1
       mark = get(Hall(ileg),'marker');
       if strcmp(mark,'none')==0
           Hmark{imark} = mark;
           imark=imark+1;
       end
   end
end

if isnan(Hleg)==0
    Hleg2 = copyobj(Hleg,hfig_new);
    hlegs=findobj(Hleg2);
    counter=1;
    imark=1;
    for ileg=1:length(hlegs)
        type = get(hlegs(ileg),'type');
        %change the fontsize for the legend text to make a bit smaller
        if strcmp(type,'text')==1
            set(hlegs(ileg),'fontsize',fontsize-2);
        end
        
        %put the markers onto the legend
        if exist('Hmark')
            Lmark = length(Hmark);
            if strcmp(type,'line')==1
                if counter==1  %mark the 2nd line
                    %set(hlegs(ileg),'marker',Hmark{imark});
                    counter=0;
                    imark=imark+1;
                else
                    counter=1;
                end
            end
        end

    end

else
    Hleg2=NaN;
end


% figure(1);
% plot(rand(100,1));grid
% title('some title'),xlabel('time');ylabel('amplitude')
% hax1=gca;
% hf2=figure(2);
% s1=subplot(211);
% pos=get(s1,'Position');
% delete(s1);
% hax2=copyobj(hax1,hf2);
% set(hax2, 'Position', pos);
   
