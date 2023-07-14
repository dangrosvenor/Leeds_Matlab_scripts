% N.B. - axis position is [left bottom width height]
Ntot = Nsub*Msub;

fontsize=14;
resizing=1;

max_tit_wrap = 30; %Max length of title before wrapping - better to avoid wrapping as changes thes
%size of the image.

savename = [fig_dir figlab];
scrsz=get(0,'ScreenSize');
%posit=[0 scrsz(4)*0.1 scrsz(3)/1.05 scrsz(4)*0.8];
%posit=[0 scrsz(4)*0.1 scrsz(3)/1.3 scrsz(4)*1.5];
posit=[0 scrsz(4) scrsz(3)*0.8 scrsz(4)];
hf=figure('position',posit,'name',figlab,'color','w');
%hf=figure('name',figlab);
%set(hf,'papertype','a4');
%hA = tight_subplot(3, 2, [.05 .00], [.1 .1], [.1 .1]);



set(hf,'paperposition',[0 0 11 8.5]);

for ifig=1:length(fig_load)
    
    [hsub_new(ifig),Hcb(ifig),pos_sub{ifig},pos_cb{ifig}] = load_fig_copy_to_subplotv2(hf,fig_dir,fig_load{ifig},ifig,fontsize,Nsub,Msub,clims{ifig});
    
    if isnan(Hcb(ifig))==1 
        no_cb=1;
    else
        no_cb=0;
     end
    
    if icbar_horiz{ifig}==1   
        if no_cb==0
            delete(Hcb(ifig));
        end
        Hcb(ifig) = colorbar('horiz');
        set(Hcb(ifig),'fontsize',fontsize);
    end

   
    
    if no_cb==0
        % -- optional --
        %Sometimes the colorbar scale gets moved to above the colorbar
        set(Hcb(ifig),'xaxislocation','bottom');
        
        
        
        
        %% adjust the sizes and colorbars
        
        
        %pos = get(hsub_new,'position');
        %axis position is [left bottom width height]
        %pos_new = [pos(1) pos(2) pos(3)*1.5 pos(4)];
        %set(hsub_new,'position',pos_new);
        
        
        %change the colorbars to be a large one underneath both. And move the plots
        %a little to reduce vertical space
        
        
        
        %% Reposition colourbars according to pos_cb_set{ifig}
        % if stretch_cb==1
        %     %the one below stretches the colorbar
        %     pos_new = [pos_cb{1}(1)+dX_cb pos_cb{1}(2)+dY_cb pos_sub{1}(3)*fcb_stretch pos_cb{1}(4)];
        % else
        %     pos_new = [pos_cb{1}(1)+dX_cb pos_cb{1}(2)+dY_cb pos_sub{1}(3) pos_cb{1}(4)];
        % end
        pos_new = [pos_cb_set{ifig}(1) pos_cb_set{ifig}(2) pos_cb_set{ifig}(3) pos_cb_set{ifig}(4)];
        set(Hcb(ifig),'position',pos_new);
        % xlabel(Hcb(1),'%');
        % % ------------------------------
        
        %axis position is [left bottom width height]
        %xlabel(Hcb(ifig),'');
        set(Hcb(ifig),'XaxisLocation','bottom');
        
        xlab = xlabel(Hcb(ifig),cbar_label{ifig},'fontsize',16);
        set(xlab,'units','normalized'); %seems to be required to get the labe lto move...
        pos = get(xlab,'position');
        %dz_tit = -0.7;
        dz_tit = 0.0;
        set(xlab,'position',[pos(1) pos(2)+dz_tit pos(3)]);
        
    end

    lb_map = lbmap(16,'brownblue');
    colormap(flipdim(lb_map,1)); 
    
    if no_cb==0 & idelete_cb{ifig}==1
        delete(Hcb(ifig));
    end
    
    h_xlab=get(hsub_new(ifig),'xlabel');
    xlab_pos = get(h_xlab,'position');
    set(h_xlab,'position',[xlab_pos(1)+dX_xlab(ifig) xlab_pos(2)+dY_xlab(ifig) xlab_pos(3)]);    
    xlab_text = get(h_xlab,'string');
    if iscell(xlab_text)==1 %for case where put a blank '' before the xlabel for spacing - seems to mess up these plots, though...
        xlabel(xlab_text{2});
    end
    set(h_xlab,'fontsize',fontsize);
    
 %% Re-size all figures by *size_fac - do this after 
 % have moved the colorbar away to make the subplot size go back to as
 % though there was no colorbar
    if resizing==1
        
        %resize the panels to remove some dead space
        %    for ifig=1:length(fig_load)
        outpos = get(hsub_new(ifig),'outerposition');
        %set(hsub_new(ifig),'outerposition',[outpos(1) outpos(2)-outpos(4)*0.1 outpos(3) outpos(4)*size_fac]);
        set(hsub_new(ifig),'outerposition',[outpos(1)+dX(ifig) outpos(2)+dY(ifig) outpos(3)*size_fac_X outpos(4)*size_fac_Y]);
        %    end
        
    end
    
    
    axes(hsub_new(ifig));
    %        xlabel('');
    %        set(gca,'ylim',[3 7]);
    %         t1 = title(textwrap({'MYD35, all r_e'},max_tit_wrap));
    t1 = title(textwrap({fig_tit{ifig}},max_tit_wrap));
    
    %Also add a title for the left and right columns
    tit_pos = get(t1,'position');
    %        str='Day';
    
    %            hlabs = annotation('textbox',[tit_pos(1) tit_pos(2)+dY_labs tit_pos(3) 1],...
    %            'Units','normalized','string', str,'linestyle','none','fontsize',fontsize);
    
    %        hcol_text = text(tit_pos(1)-0.04,tit_pos(2)+dY_labs*4,str);
    %        set(hcol_text,'fontsize',18);
    
    %xticks = get(hsub_new(ifig),'XTick');
    %xticks = xticks * xscale;
    %h_xlab = get(hsub_new(ifig),'Xlabel');
    %hlab_text = get(h_xlab,'string');
    %hlab_text{2} = 
    %xticks = get(hsub_new(ifig),'XTick');
    
    if length(strmatch(ioverride_xticks{ifig},''))==0
        set(hsub_new(ifig),'XTick',ioverride_xticks{ifig});
    end
    if length(strmatch(ioverride_yticks{ifig},''))==0
        set(hsub_new(ifig),'YTick',ioverride_yticks{ifig});
    end
    if idelete_xlab{ifig}==1;
        xlabel('');
    end
    
    if idelete_ylab{ifig}==1;
        ylabel('');
    end
    if idelete_xtick_labs{ifig}==1
        set(hsub_new(ifig),'XTickLabel','');
    end
    
   

    
    % -------------------------------------
    
    %% label with a,b,c, etc.
    if Ntot>1
        abc = {'abcdefghijklmnopqrstuvwxyz'};
        %        axes(hsub_new(ifig)); %make the required plot current
        outpos = get(hsub_new(ifig),'position');
        %[left bottom width height]
        
        %text is written relative to current plot
        str = ['(' abc{1}(ifig) ')'];
        %        text(0,outpos(4)+dY,['(' abc{1}(ifig) ')']);
        
        %but use this instead for relative to the figure - PROBLEM - this
        %uses Java and produces poor quality text resolution upon printing.
        %Apparently there is no solution to this...
        %Annotation below works better.
        %         hlabs = uicontrol('Style', 'text',...
        %        'String', str ,... %replace something with the text you want
        %        'Units','normalized',...
        %        'Position', [outpos(1)-0.02 outpos(2)+outpos(4)-0.1 0.03 0.03],...
        %        'backgroundcolor',[1 1 1],'Fontname','Helvetica');
        
        
        
        %[left bottom width height]
        hlabs = annotation('textbox',[min(max(outpos(1)+dX_labs,0),1) max(min(outpos(2)+outpos(4)+dY_labs,1),0) 0.03 0.03],...
            'Units','normalized','string', str,'linestyle','none','fontsize',fontsize);
        
        
        %         if ifig==1
        %             halfwidth=outpos(3)/2;
        %             str = 'Day';
        %             hlabs = annotation('textbox',[min(max(outpos(1),0),1) max(min(outpos(2)+outpos(4)+dY_labs*3,1),0) 0.03 0.03],...
        %             'Units','normalized','string', str,'linestyle','none','fontsize',fontsize);
        %         end
        
        
        
        
        
    end
    
    
    
    
end

%set(gcf,'renderer','zbuffer'); %painters is the one that produces scalable .eps
%However, it is fairly slow and memory hungry. Plus the vectorized pcolor
%style maps have the anti-aliasing lines - here are forcing to produce a
%high resolution bitmap (but still in .eps format)

datestr_now = saveas_ps_fig_emf(gcf,savename,'',0,0,0);




