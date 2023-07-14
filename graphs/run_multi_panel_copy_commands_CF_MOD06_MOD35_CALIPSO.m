savename = [fig_dir figlab];
%scrsz=get(0,'ScreenSize');
%posit=[0 scrsz(4)*0.1 scrsz(3)/1.05 scrsz(4)*0.8];
%posit=[0 scrsz(4)*0.1 scrsz(3)/1.3 scrsz(4)*1.5];
%hf=figure('position',posit,'name',figlab);
hf=figure('name',figlab);
%set(hf,'papertype','a4');
%hA = tight_subplot(3, 2, [.05 .00], [.1 .1], [.1 .1]);

fontsize=16;

set(hf,'paperposition',[0 0 8.5 11]);




% ifig=ifig+1;
% fig_load = 'Low_Cloud_Fraction_COSP-CALIPSO_ANNUAL,_2006-2010_DAYTIME_for_2007_to_2010__mean=0.3_20130718T185000.fig';
% hf_load = open([fig_dir fig_load]);
% title(''); %remove the title
% hs = subplot_spaceplots(3,2,ifig,'parent',hf);
% %hs = hA(1);
% [hsub_new(ifig),Hcb,pos_sub,pos_cb]=fig2subplot(hf_load,hs);
% temp = increase_font_size_map_figures_func(hsub_new(ifig),16,0);
% close(gcf);




%ifig=ifig+1;
%fig_load = 'Low_Cloud_Fraction_COSP-CALIPSO_ANNUAL,_2006-2010_NIGHTTIME_for_2007_to_2010__mean=0_20130718T184433.fig';

for ifig=1:length(fig_load)

[hsub_new(ifig),Hcb(ifig),pos_sub{ifig},pos_cb{ifig}] = load_fig_copy_to_subplot(hf,fig_dir,fig_load{ifig},ifig,fontsize);
%hf_load = open([fig_dir fig_load]);
%title(''); %remove the title
%hs = subplot_spaceplots(3,2,ifig,'parent',hf);
%hs = hA(2);
%[hsub_new(ifig),Hcb2,pos_sub2,pos_cb2]=fig2subplot(hf_load,hs);
%temp = increase_font_size_map_figures_func(hsub_new(ifig),16,0);
%close(gcf);

end
% 
% ifig=ifig+1;
% fig_load = 'Low_Cloud_Fraction_COSP-CALIPSO_ANNUAL,_2006-2010_NIGHTTIME_for_2007_to_2010__mean=0_20130718T184433.fig';
% hf_load = open([fig_dir fig_load]);
% title(''); %remove the title
% hs = subplot_spaceplots(3,2,ifig,'parent',hf);
% %hs = hA(2);
% [hsub_new(ifig),Hcb2,pos_sub2,pos_cb2]=fig2subplot(hf_load,hs);
% temp = increase_font_size_map_figures_func(hsub_new(ifig),16,0);
% close(gcf);
% 
% ifig=ifig+1;
% fig_load = 'Low_Cloud_Fraction_COSP-CALIPSO_ANNUAL,_2006-2010_NIGHTTIME_for_2007_to_2010__mean=0_20130718T184433.fig';
% hf_load = open([fig_dir fig_load]);
% title(''); %remove the title
% hs = subplot_spaceplots(3,2,ifig,'parent',hf);
% %hs = hA(2);
% [hsub_new(ifig),Hcb2,pos_sub2,pos_cb2]=fig2subplot(hf_load,hs);
% temp = increase_font_size_map_figures_func(hsub_new(ifig),16,0);
% close(gcf);
% 
% ifig=ifig+1;
% fig_load = 'Low_Cloud_Fraction_COSP-CALIPSO_ANNUAL,_2006-2010_NIGHTTIME_for_2007_to_2010__mean=0_20130718T184433.fig';
% hf_load = open([fig_dir fig_load]);
% title(''); %remove the title
% hs = subplot_spaceplots(3,2,ifig,'parent',hf);
% %hs = hA(2);
% [hsub_new(ifig),Hcb2,pos_sub2,pos_cb2]=fig2subplot(hf_load,hs);
% temp = increase_font_size_map_figures_func(hsub_new(ifig),16,0);
% close(gcf);
% 
% ifig=ifig+1;
% fig_load = 'Low_Cloud_Fraction_COSP-CALIPSO_ANNUAL,_2006-2010_NIGHTTIME_for_2007_to_2010__mean=0_20130718T184433.fig';
% hf_load = open([fig_dir fig_load]);
% title(''); %remove the title
% hs = subplot_spaceplots(3,2,ifig,'parent',hf);
% %hs = hA(2);
% [hsub_new(ifig),Hcb2,pos_sub2,pos_cb2]=fig2subplot(hf_load,hs);
% temp = increase_font_size_map_figures_func(hsub_new(ifig),16,0);
% close(gcf);

%% adjust the sizes and colorbars

%resize the panels to remove some dead space
size_fac = 1.3;
for ifig2=1:length(fig_load)
    outpos = get(hsub_new(ifig2),'outerposition');
    %set(hsub_new(ifig2),'outerposition',[outpos(1) outpos(2)-outpos(4)*0.1 outpos(3) outpos(4)*size_fac]);
    set(hsub_new(ifig2),'outerposition',[outpos(1) outpos(2) outpos(3) outpos(4)*size_fac]);
end



%pos = get(hsub_new,'position');
%axis position is [left bottom width height]
%pos_new = [pos(1) pos(2) pos(3)*1.5 pos(4)];
%set(hsub_new,'position',pos_new);


%change the colorbars to be a large one underneath both. And move the plots
%a little to reduce vertical space

dY = 0.02; %relative to papersize
dY_cb = 0.05; %distance to move the colorbars upwards

%top two plots
delete(Hcb(2)); %delete the right colorbar
%axis position is [left bottom width height]
pos_new = [pos_cb{1}(1) pos_cb{1}(2)+dY_cb pos_sub{1}(3)*2.2 pos_cb{1}(4)];
set(Hcb(1),'position',pos_new);
%xlabel(Hcb(1),'Low Cloud Fraction');
for ifig2=1:2
    outpos = get(hsub_new(ifig2),'outerposition');
    set(hsub_new(ifig2),'outerposition',[outpos(1) outpos(2)-dY outpos(3) outpos(4)]);
end

%middle two plots
delete(Hcb(4)); %delete the right colorbar
%axis position is [left bottom width height]
pos_new = [pos_cb{3}(1) pos_cb{3}(2)+dY_cb pos_sub{3}(3)*2.2 pos_cb{3}(4)];
set(Hcb(3),'position',pos_new);
%xlabel(Hcb(3),'Middle Cloud Fraction');
for ifig2=3:4
    outpos = get(hsub_new(ifig2),'outerposition');
    set(hsub_new(ifig2),'outerposition',[outpos(1) outpos(2)-dY outpos(3) outpos(4)]);
end

%bottom two plots
delete(Hcb(6)); %delete the right colorbar
%axis position is [left bottom width height]
pos_new = [pos_cb{5}(1) pos_cb{5}(2)+dY_cb pos_sub{5}(3)*2.2 pos_cb{5}(4)];
set(Hcb(5),'position',pos_new);
%xlabel(Hcb(5),'High Cloud Fraction');
for ifig2=5:6
    outpos = get(hsub_new(ifig2),'outerposition');
    set(hsub_new(ifig2),'outerposition',[outpos(1) outpos(2)-dY outpos(3) outpos(4)]);
end

%htit(1) = title(hsub_new(1),'Day');
%htit(2) = title(hsub_new(2),'Night');

%spaceplots(hf,[0 0 0 0], [.02 .02]);

%% label with a,b,c, etc.
abc = {'abcdefghijklmnopqrstuvwxyz'};
for ifig2=1:length(fig_load)
%        axes(hsub_new(ifig2)); %make the required plot current
        outpos = get(hsub_new(ifig2),'outerposition');
        %text is written relative to current plot
        str = ['(' abc{1}(ifig2) ')'];
%        text(0,outpos(4)+dY,['(' abc{1}(ifig2) ')']);
        
        %but use this instead for relative to the figure - PROBLEM - this
        %uses Java and produces poor quality text resolution upon printing.
        %Apparently there is no solution to this...
        %Annotation below works better.        
%         hlabs = uicontrol('Style', 'text',...
%        'String', str ,... %replace something with the text you want
%        'Units','normalized',...
%        'Position', [outpos(1)-0.02 outpos(2)+outpos(4)-0.1 0.03 0.03],...
%        'backgroundcolor',[1 1 1],'Fontname','Helvetica'); 
   
   hlabs = annotation('textbox',[outpos(1)-0.02 outpos(2)+outpos(4)-0.1 0.03 0.03],...
       'Units','normalized','string', str,'linestyle','none','fontsize',fontsize);
   
   
        
end



set(gcf,'renderer','zbuffer'); %painters is the one that produces scalable .eps
%However, it is fairly slow and memory hungry. Plus the vectorized pcolor
%style maps have the anti-aliasing lines - here are forcing to produce a
%high resolution bitmap (but still in .eps format)
datestr_now = saveas_ps_fig_emf(gcf,savename,'',0,1,0);





