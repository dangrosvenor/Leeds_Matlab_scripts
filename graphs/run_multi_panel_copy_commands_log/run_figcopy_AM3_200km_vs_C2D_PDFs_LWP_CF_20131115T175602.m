%N.B. Need to set color limits before this otherwise colobar will not
%change - specify in the script with the figure names in
Nsub = 1;
Msub = 2;  %No. of panels Nsub = no. vert, Msub = no. horiz

Ntot = Nsub*Msub;

fontsize=16;
resizing=0;
stretch_cb=0; %whether to stretch the CB across the bottom of the plot

%positive is up the page, negative down
fcb_stretch = 1.7; %2.2
size_fac = 0.27;  %0.78; %0.85;
dY = 0.25; %0.105;  %0.09; %0.02; %distance to move the subplot upwards relative to papersize
dX = -0.10;
even_only_dX = 1; %flag to tell it to only move the odd numbered plot by dX
dY_cb = 0.36; %0.125;  %0.11; %distance to move the colorbars upwards
dY_spacing = 0.05; %vertical spacing to add between plots

%move the a,b,c labels in x and y by:-
dX_labs = -0.02;
dX_labs = 0.02;
dX_labs = 0; %-0.1;

dY_labs = -0.1;
dY_labs = -0.05;
dY_labs = -0.05; %0.02; 

savename = [fig_dir figlab];
%scrsz=get(0,'ScreenSize');
%posit=[0 scrsz(4)*0.1 scrsz(3)/1.05 scrsz(4)*0.8];
%posit=[0 scrsz(4)*0.1 scrsz(3)/1.3 scrsz(4)*1.5];
%hf=figure('position',posit,'name',figlab);
hf=figure('name',figlab);
%set(hf,'papertype','a4');
%hA = tight_subplot(3, 2, [.05 .00], [.1 .1], [.1 .1]);



set(hf,'paperposition',[0 0 11 8.5]);

for ifig=1:length(fig_load)

[hsub_new(ifig),Hcb(ifig),pos_sub{ifig},pos_cb{ifig}] = load_fig_copy_to_subplotv2(hf,fig_dir,fig_load{ifig},ifig,fontsize,Nsub,Msub,clims{ifig});

% -- optional --
%Sometimes the colorbar scale gets moved to above the colorbar
%set(Hcb(ifig),'xaxislocation','bottom');

end

%% adjust the sizes and colorbars


if resizing==1

    %resize the panels to remove some dead space
    for ifig2=1:length(fig_load)
        outpos = get(hsub_new(ifig2),'outerposition');
        %set(hsub_new(ifig2),'outerposition',[outpos(1) outpos(2)-outpos(4)*0.1 outpos(3) outpos(4)*size_fac]);
        set(hsub_new(ifig2),'outerposition',[outpos(1) outpos(2) outpos(3) outpos(4)*size_fac]);
    end

end

%pos = get(hsub_new,'position');
%axis position is [left bottom width height]
%pos_new = [pos(1) pos(2) pos(3)*1.5 pos(4)];
%set(hsub_new,'position',pos_new);


%change the colorbars to be a large one underneath both. And move the plots
%a little to reduce vertical space



%% Delete one of the colorbars and make the other one go across the whole
%% page below the top two subplots

% % ---- optional settings ------
delete(Hcb(2)); %delete the right colorbar
%axis position is [left bottom width height]


if stretch_cb==1
    %the one below stretches the colorbar
    pos_new = [pos_cb{1}(1) pos_cb{1}(2)+dY_cb pos_sub{1}(3)*fcb_stretch pos_cb{1}(4)];
else
    pos_new = [pos_cb{1}(1) pos_cb{1}(2)+dY_cb pos_sub{1}(3) pos_cb{1}(4)];
end
 set(Hcb(1),'position',pos_new);
 xlabel(Hcb(1),'Frequency','fontsize',fontsize-4);
% % ------------------------------



    for ifig2=1
        if resizing==1
            dX_val = dX;
            if rem(ifig2,2)~=0 %if an odd number
                if even_only_dX ==1
                    dX_val = 0;
                end
            end
            outpos = get(hsub_new(ifig2),'outerposition');
            set(hsub_new(ifig2),'outerposition',[outpos(1)+dX_val outpos(2)+dY outpos(3) outpos(4)]);
        end
        axes(hsub_new(ifig2));
%        xlabel('');
%        set(gca,'ylim',[3 7]);
%         title('Day');
    end


if Ntot >= 2
     for ifig2=2
        if resizing==1
           dX_val = dX;
            if even_only_dX ==1
                if rem(ifig2,2)~=0 %if an odd number    
                    dX_val = 0;
                end
            end   

            outpos = get(hsub_new(ifig2),'outerposition');
            set(hsub_new(ifig2),'outerposition',[outpos(1)+dX_val outpos(2)+dY outpos(3) outpos(4)]);
        end
        axes(hsub_new(ifig2));
        
% ---- optional settings ------        
%        xlabel('');
%        ylabel('$\theta_{1 km}$','interpreter','latex');
%        set(gca,'ylim',[225 300]);
%         title('Night');
     end
     

     


end
    

if Ntot>=3
    dY = dY + dY_spacing;  %add on dY_spacing before odd numbered plots
    dY_cb = dY_cb + dY_spacing;  %add on dY_spacing before odd numbered plots
    
    % ---- optional settings ------
    delete(Hcb(4)); %delete the right colorbar
    %axis position is [left bottom width height]
    if stretch_cb==1
        pos_new = [pos_cb{3}(1) pos_cb{3}(2)+dY_cb pos_sub{3}(3)*fcb_stretch pos_cb{3}(4)];
    else
        pos_new = [pos_cb{3}(1) pos_cb{3}(2)+dY_cb pos_sub{3}(3) pos_cb{3}(4)];
    end
    set(Hcb(3),'position',pos_new);
    xlabel(Hcb(3),'COSP CALIPSO CF bias','fontsize',fontsize-4);
% ------------------------------


    for ifig2=3
        if resizing==1
           dX_val = dX;
            if even_only_dX ==1
                if rem(ifig2,2)~=0 %if an even number     
                    dX_val = 0;
                end
            end            
            outpos = get(hsub_new(ifig2),'outerposition');
            set(hsub_new(ifig2),'outerposition',[outpos(1)+dX_val outpos(2)+dY outpos(3) outpos(4)]);
        end
        axes(hsub_new(ifig2));
        
% ---- optional settings ------
%        xlabel('');
%        ylabel('$\theta_{2 km}$','interpreter','latex');
%        set(gca,'ylim',[225 300]);
%         title('MOD06, no height screening');
          title(' '); %Add a space for the title, otherwise figure sizes are different...
% ------------------------------         
    end
    
end
    

if Ntot>=4
     for ifig2=4
        if resizing==1
            dX_val = dX;
            if even_only_dX ==1
                if rem(ifig2,2)~=0 %if an even number   
                    dX_val = 0;
                end
            end            
            outpos = get(hsub_new(ifig2),'outerposition');
            set(hsub_new(ifig2),'outerposition',[outpos(1)+dX_val outpos(2)+dY outpos(3) outpos(4)]);
        end
        axes(hsub_new(ifig2));
        
% ---- optional settings ------        
 %       xlabel('');
  %      set(gca,'ylim',[0.014 0.019]);
%         title('MOD06, CTH < 3.25km');
          title(' '); %Add a space for the title, otherwise figure sizes are different...-
    end






%htit(1) = title(hsub_new(1),'Day');
%htit(2) = title(hsub_new(2),'Night');

%spaceplots(hf,[0 0 0 0], [.02 .02]);

end

if Ntot>=5
    dY = dY + dY_spacing;  %add on dY_spacing before odd numbered plots
     dY_cb = dY_cb + dY_spacing;  %add on dY_spacing before odd numbered plots
     
    % ---- optional settings ------
    delete(Hcb(6)); %delete the right colorbar
    %axis position is [left bottom width height]
    
    stretch_cb=1;
    if stretch_cb==1
        pos_new = [pos_cb{5}(1) pos_cb{5}(2)+dY_cb pos_sub{5}(3)*fcb_stretch pos_cb{5}(4)];
    else
        pos_new = [pos_cb{5}(1) pos_cb{5}(2)+dY_cb pos_sub{5}(3) pos_cb{5}(4)];
    end
    set(Hcb(5),'position',pos_new);
    xlabel(Hcb(5),'CAMCLUBBv2','fontsize',fontsize-4);
% ------------------------------


 for ifig2=5
        if resizing==1
            dX_val = dX;
            if even_only_dX ==1
                if rem(ifig2,2)~=0 %if an even number   
                    dX_val = 0;
                end
            end              
            outpos = get(hsub_new(ifig2),'outerposition');
            set(hsub_new(ifig2),'outerposition',[outpos(1)+dX_val outpos(2)+dY outpos(3) outpos(4)]);
        end
        axes(hsub_new(ifig2));
        
% ---- optional settings ------
%        xlabel('');
%        ylabel('$\theta_{2 km}$','interpreter','latex');
%        set(gca,'ylim',[225 300]);
%         title('MOD35 - CALIPSO, CTH < 3.25km');
          title(' '); %Add a space for the title, otherwise figure sizes are different...
% ------------------------------         
 end
 
end


if Ntot>=6
    
   

  for ifig2=6
      
        % ---- optional settings ------
%    delete(Hcb(6)); %delete the right colorbar
    %axis position is [left bottom width height]
    
%     stretch_cb=0;
%     if stretch_cb==1
%         pos_new = [pos_cb{ifig2}(1) pos_cb{ifig2}(2)+dY_cb pos_sub{ifig2}(3)*fcb_stretch pos_cb{ifig2}(4)];
%     else
%         pos_new = [pos_cb{ifig2}(1) pos_cb{ifig2}(2)+dY_cb pos_sub{ifig2}(3) pos_cb{ifig2}(4)];
%     end
%     set(Hcb(ifig2),'position',pos_new);
%    xlabel(Hcb(ifig2),'%');
% ------------------------------


        if resizing==1
            dX_val = dX;
            if even_only_dX ==1
                if rem(ifig2,2)~=0 %if an even number   
                    dX_val = 0;
                end
            end              
            outpos = get(hsub_new(ifig2),'outerposition');
            set(hsub_new(ifig2),'outerposition',[outpos(1)+dX_val outpos(2)+dY outpos(3) outpos(4)]);
        end
        axes(hsub_new(ifig2));
        
% ---- optional settings ------
%        xlabel('');
%        ylabel('$\theta_{2 km}$','interpreter','latex');
%        set(gca,'ylim',[225 300]);
%         title('MOD06 - CALIPSO, CTH < 3.25km');
          title(' '); %Add a space for the title, otherwise figure sizes are different...
% ------------------------------         
  end
    
end



% -------------------------------------

%% label with a,b,c, etc.
if Ntot>1
    abc = {'abcdefghijklmnopqrstuvwxyz'};
    for ifig2=1:length(fig_load)
        %        axes(hsub_new(ifig2)); %make the required plot current
        outpos = get(hsub_new(ifig2),'position');
        %[left bottom width height]

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



        %[left bottom width height]
        hlabs = annotation('textbox',[max(outpos(1)+dX_labs,0) min(outpos(2)+outpos(4)+dY_labs,1) 0.03 0.03],...
            'Units','normalized','string', str,'linestyle','none','fontsize',fontsize);



    end

end


%set(gcf,'renderer','zbuffer'); %painters is the one that produces scalable .eps
%However, it is fairly slow and memory hungry. Plus the vectorized pcolor
%style maps have the anti-aliasing lines - here are forcing to produce a
%high resolution bitmap (but still in .eps format)

datestr_now = saveas_ps_fig_emf(gcf,savename,'',0,1,0);





