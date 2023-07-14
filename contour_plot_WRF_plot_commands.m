





maxCov=max(maxC);
minCov=min(minC);

if isnan(maxCov); maxCov=0; end
if isnan(minCov); minCov=0; end





%need to make sure highest contour is higher than highest data value for colorbarf

if i2d==1
    timesTH(1).t=Grid.Y1/1000;
    xlabelstr='Horizontal Distance (km)';
elseif i2d==2
    xlabelstr='Horizontal Distance (km)';
elseif i2d==0
    if iutc==1
        xlabelstr='UTC Time (hrs)';
        xlabelstr='Time (hrs)';

    else
        xlabelstr='Local Time (hrs)';
    end
end


if isamescale==1
    iminovr=1;
    imaxovr=1;
    if dlogflag==1
        mincovOvr = dlog(minVal,dlogmin);
        maxcovOvr = dlog(maxVal,dlogmin);
    elseif logflag==1
        mincovOvr = log10(minVal);
        maxcovOvr = log10(maxVal);
    else
        mincovOvr = minVal;
        maxcovOvr = maxVal;
    end
end

for i=1:nplots2d

    if notsame==1
        maxCov=maxC(i);
        minCov=minC(i);
    end

    if iminovr(i)==1    %imincovovr==
        minCov=mincovOvr(i);
    else
        if abs(minC(i))>=100 | abs(maxC(i))>=100 %if length(num2str(round(maxC(i)))) - length(num2str(round(minC(i)))) == 0
            nsigfig = 3
        elseif abs(minC(i))>=10 | abs(maxC(i))>=10
            nsigfig = 2;
        else
            nsigfig = 1;
        end
        mincovOvr(i) = sigfig(minC(i),nsigfig);   %Dan - added to make nicer numbers on contours and colorbar
        minCov=mincovOvr(i);
        iminovr(i)=1;
    end

    if imaxovr(i)==1    %imaxcovovr==
        maxCov=maxcovOvr(i);
    else
        maxcovOvr(i) = sigfig(maxC(i),nsigfig);
        maxCov=maxcovOvr(i);
        imaxovr(i)=1;
    end

    minc=minCov-0.0*abs(minCov);
    maxc=maxCov+abs(maxCov*0.0);

    if minc>maxc; m=minc; minc=maxc; maxc=m; end

    %     fixmin=fix2(minc,abs(round(log10(min(abs(minc))))));
    %     dd=(maxc-minc)/ncont;
    %     dfix=fix2(dd,ceil(abs(log10(min(abs(dd))))));
    %     conts=[fixmin:dfix:maxc];
    %conts=[minc:(maxc-minc)/ncont:maxc];
    %conts=round2(conts,abs(round(log10(min(abs(conts)))))+1);

    if logflag==0 & dlogflag==0
        dd=(maxc-minc)/ncont;
        %dfix=fix2(dd,ceil(abs(log10(min(abs(dd))))));

        dfix=sigfig(dd,0);
        if dfix==0;
            dfix=sigfig(dd,1);
        end
        if minc>=100
            fixmin=sigfig(minc,sig);
        else
            fixmin=sigfig(minc,sig);
        end
        if abs(fixmin)<dfix/100; fixmin=0; end
        conts=[fixmin:dfix:maxc];
        iszero=zeros([1 length(conts)+1]); %flag to say whether value is exactly zero

        if length(conts)==0;
            conts(1)=0;
            conts(2)=0;
        end

        if sign(conts(1))~=sign(conts(end))
            [zeromin izero]=min(abs(conts));
            if abs(conts(izero))<dfix/2.2
                conts(izero)=0;
                iszero(izero)=1; %flag to say that value is exactly zero
            else
                if conts(izero)<0
                    conts(izero+2:end+1)=conts(izero+1:end);
                    conts(izero+1)=0;
                    iszero(izero+1)=1;
                else
                    conts(izero+1:end+1)=conts(izero:end);
                    conts(izero)=0;
                    iszero(izero)=1;
                end
            end
        end
    elseif dlogflag==1
        conts=[minc:(maxc-minc)/ncont:maxc];
        iszero=zeros([1 length(conts)+1]);

        unlog=idlog(conts,dlogmin);
        if sign(unlog(1))~=sign(unlog(end))
            [zeromin izero]=min(abs(unlog));
            if abs(unlog(izero))<((maxc-minc)/ncont)/2.2
                unlog(izero)=0;
                iszero(izero)=1;
            else
                if unlog(izero)<0
                    unlog(izero+2:end+1)=unlog(izero+1:end);
                    unlog(izero+1)=0;
                    iszero(izero+1)=1;
                else
                    unlog(izero+1:end+1)=unlog(izero:end);
                    unlog(izero)=0;
                    iszero(izero)=1;
                end
            end
        end

        unlog=sigfig(unlog,sig);
        conts=dlog(unlog,dlogmin);

    else
        conts=[minc:(maxc-minc)/ncont:maxc]; %logflag==1
        iszero=zeros([1 length(conts)+1]);
    end


    if length(conts)==0
        conts=[0 1];
    end

    %     ac=find(conts>0);
    %     ac2=find(conts<0);
    %     if length(ac)>0 & ac(1)>1
    %         ac=ac(1);
    %         conts=[conts(1:ac-1) 0 conts(ac:end)];
    %     end

    if iovride_conts==1
        conts=conts_ovr;
    end



    % put in the extra bits at the side to get the filled contour/filled colorbar right
    psame=repmat(pdat(i).p(:,pend(i)),[1 npend(i)-pend(i)-nend(i)-1]);
    %    psame(isnan(psame))=0;
    pdat(i).p(:,pend(i)+1:npend(i)-nend(i)-1)=psame;


    mindat=minCov*(1-0.01*sign(minCov));
    maxdat=maxCov*(1+0.01*sign(maxCov));

    switch right_side_extra_bits
        case 1
            pdat(i).p(1:half(i),npend(i)-nend(i):npend(i))=mindat;
            pdat(i).p(half(i)+1:end,npend(i)-nend(i):npend(i))=maxdat;
        case 2
            pdat_dummy = pdat(i).p;
            pdat_dummy(1:half(i),npend(i)-nend(i):npend(i))=mindat;
            pdat_dummy(half(i)+1:end,npend(i)-nend(i):npend(i))=maxdat;
    end




    if subplotting==0
        h(isub).h=subplot(a,b,isub);
    end


    if izovr==0;
        zz(i).z=z(izmin:izmax);
    end

    %pcolor(timesTH,0.62+z(izmin:izmax)./1000,pdat);
    if ilem==1;
        height=add_ground_height+zz(i).z./1000;
    else
        height=zz(i).z;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%    draw the plot    %%%%%%%%%%%%%%%%%%%%%%%%%

    if icont==1
        if clines==1
            %                [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts); %old version
            %                 [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts,'linestyle','none');
            %                 [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts,'linestyle','none');
            %                 hold on

            [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts,'linestyle','none');
            hold on
            %[cbfA(i).cpos cbfB(i).cpos]=contourf(timesTH(i).t,height,pdat(i).p,conts(conts<0),'k','linestyle',':');
            %[cbfA(i).cneg cbfB(i).cneg]=contourf(timesTH(i).t,height,pdat(i).p,conts(conts>=0),'k','linestyle','-');
            %changed the above two contourf calls to contour calls 26/03/09
            [cbfA(i).cpos cbfB(i).cpos]=contour(timesTH(i).t,height,pdat(i).p,conts(conts<0),'k','linestyle',':');
            [cbfA(i).cneg cbfB(i).cneg]=contour(timesTH(i).t,height,pdat(i).p,conts(conts>=0),'k','linestyle','-');
            icont_neg=1;
            %hold on
            % [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts,'linestyle','none');



        else
            [cbfA(i).c cbfB(i).c]=contourf(timesTH(i).t,height,pdat(i).p,conts,'linestyle','none');
            cbfA(i).cpos = cbfA(i).c;
            cbfB(i).cpos = cbfB(i).c;
        end
    else
        [cbfA(i).c]=pcolor(timesTH(i).t,height,pdat(i).p);shading interp;
    end




    if icont_extra==1
        hold on
        [c2,h2]=contour(timesTH(i).t(1:length(xinds)),add_ground_height+zz(i).z/1000 ...
            ,f/1000*sum(TwoDDan(idir).Q(izmin:izmax,xinds,4:6),3),[0:5:20],'w');
        set(h2,'linewidth',1.5);
        clabel(c2,h2,'color','w','fontsize',max([fsize-8 6]));
        %        clabel(c2,h2,'labelspacing',72,'color','w');
    end



    if isquare==1
        axis square
    end


    if icont==0 %orginally this
        %    if dan_test==1 %%% dan test %%%
        cax=get(gca,'clim');
        if iminovr==1
            cax(1)=mincovOvr;
        end
        if imaxovr==1
            cax(2)=maxcovOvr;
        end
        caxis(cax);
        %        hc=colorbar;
        %        set(hc,'fontsize',fsize-3);

    end

    %hc=colorbarf(cbfA,cbfB);



    %xti=set(h(i).h,'xticklabels',timestxt);

    set(h(isub).h,'fontsize',fsize);

    if onexlabel==0 | (isub==nplots2d & onexlabel==1 & subplotting==0) | (isub==length(idirs) & onexlabel==1 & subplotting==1)
        xlabel(xlabelstr);
    end

    if izovr~=2
        ylabel('Height (km)');
    else
        ylabel(ylabelstr);
    end

    if subplotting==1 & isub==1
        title(tit(i).tit,'fontsize',fsize-4);
    elseif subplotting==0
        title(tit(i).tit,'fontsize',fsize-4);
    end



    if length(cbfA(i).c)==0
        nocbar=1;
        normcbar=1;
    end

    if icont==0
        normcbar=1;
    end
    %caxis(h(i).h,[minCov maxCov*1.05]);





    if exist('idirs');
        if isub==length(idirs)
            %    if isamescale==1 &  isub==length(idirs)

            pos1=get(h(1).h,'position'); %[left bottom width height]
            posend=get(h(end).h,'position');

            height_pos=pos1(2)-posend(2)+pos1(4);
            pos=[posend(1) posend(2) posend(3)+0.15 height_pos];
        else
            pos=[0 0 1 1];
        end
    end



    if isamescale==1 &  isub==length(idirs)
        if normcbar==1 & bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
            hc=colorbar( 'peer' , h(isub).h );
        elseif normcbar==1 & bigcbar==0 | (bigcbar==1 & isub==length(idirs)) | lememm==1
            axdan=axes('position',pos,'visible','off');
            %colbar=colorbar; %if colorbar already in place then colorbarf will replace it
            %set(colbar,'tag','Colorbar');
            hc=colorbarf(cbfA(i).c,cbfB(i).c);


        end
    else
        if nocbar~=1 & icont==1
            if normcbar==0 & bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
                %%%%%   usual one %%%%%%%%
                hc=colorbarf_30Sep09_2(cbfA(i).c,cbfB(i).c); %NOTE setting colourbar settings to be same as second graph so make sure are the same!!!!!!!!!!!!!!
                %                hc=colorbarf(cbfA(i).c,cbfB(i).c); %NOTE setting colourbar settings to be same as second graph so make sure are the same!!!!!!!!!!!!!!
                %    hc=cbarf(conts,conts,'vertical','linear');
            elseif normcbar==0 &(bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
                axdan=axes('position',pos,'visible','off')
                %   colbar=colorbar;
                %   set(colbar,'tag','Colorbar');
                hc=colorbarf(cbfA(i).c,cbfB(i).c);
                %

            end
        end

        if normcbar==1 & bigcbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
            hc=colorbar( 'peer' , h(isub).h );
        end

    end

    %   if normcbar==0

    if (logflag==1 | dlogflag==1 ) & nocbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
        %re-label colorbar ticks if have log settings
        clear ctickstr

        ctick=get(hc,'yticklabel');
        if normcbar==1 %since for normal cbar ctick returned as matrix of characters
            ctick2=ctick;
            clear ctick
            for j=1:size(ctick2,1)
                ctick{j}=ctick2(j,:);
            end
            jinds=1:length(ctick);
        else
            jinds=2:length(ctick)-1;
        end

        ctickstr(1,1)=' ';
        for j=jinds
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            nu=str2num(ctick{j});



            %te=num2str(10^nu,'%2.2f');       %'%2.2e');
            if dlogflag==0
                if 10^nu>99
                    te=num2str(sigfig(10^nu - offset ,3));
                else
                    te=num2str(sigfig(10^nu - offset,2));
                end
            else
                %                     if idlog(nu,dlogmin)>99
                %                         te=num2str(sigfig(idlog(nu,dlogmin) - offset ,4));
                %                     else
                %                         te=num2str(sigfig(idlog(nu,dlogmin) - offset,4));
                %                     end
                te=num2str(sigfig(idlog(nu,dlogmin) - offset,2),4);
                if normcbar==1
                    if iszero(j-1)==1  %if is suppossed to be exactly zero then set to zero
                        te='0.0';
                    end
                else
                    if iszero(j-1)==1  %if is suppossed to be exactly zero then set to zero
                        te='0.0';
                    end
                end

            end

            ctickstr(j,1:length(te))=te;
        end

        set(hc,'yticklabel',ctickstr);

        %         add=str2num(ctick{end-1})/50;
        %
        %         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
        %
        %         for i=2:length(ctick)-1
        %             cticknums(i)=str2num(ctick{i});
        %         end
        %         text(  ones( length(cticknums),1 )*1.05,cticknums+add,ctickstr, 'fontsize',fsize-6  );

    elseif (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) ) %if (logflag==1 | dlogflag==1 ) & nocbar=...
        %re-label colorbar if not log plot
        clear ctickstr
        ctick=get(hc,'yticklabel');

        if normcbar==1 %since for normal cbar ctick returned as matrix of characters
            ctick2=ctick;
            clear ctick
            for j=2:size(ctick2,1)+1
                ctick{j}=ctick2(j-1,:);
            end
            jinds=2:length(ctick);
        else
            jinds=2:length(ctick)-1;
        end


        ctick2=ctick;
        clear ctick
        for j=jinds
            ctick(j-1)=str2num(ctick2{j});
        end
        ctickstr(1,:)=' ';

        for j=2:length(ctick)+1
            %te=strcat('10^','{',num2str(ctick(j)),'}' );
            %te=num2str(ctick(j),'%2.2e');
            te=num2str(sigfig(ctick(j-1),sig));
            if  normcbar==0
                ctickstr(j,1:length(te))=te;
            else
                ctickstr(j-1,1:length(te))=te;
            end
        end

        set(hc,'yticklabel',ctickstr);

        %         add=ctick(end)/50;
        %
        %         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
        %         text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize  );





        %set(hc,'fontsize',fsize-6);



        %else



        %     if logflag==1 & nocbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
        %         clear ctickstr
        %         ctick=get(hc,'ytick');
        %         for j=1:length(ctick)
        %             %te=strcat('10^','{',num2str(ctick(j)),'}' );
        %             %te=num2str(10^ctick(j),'%2.2g');
        %             if nu>99
        %                 te=num2str(sigfig(10^ctick(j) - offset,2));
        %             else
        %                 te=num2str(sigfig(10^ctick(j) - offset,2));
        %             end
        %
        %             ctickstr(j,1:length(te))=te;
        %         end
        %
        %         set(hc,'yticklabel','');
        %
        %         add=ctick(end)/50;
        %
        %         set(hf,'currentaxes',hc); %this also allows you to use xlabel, ylabel and title for colorbar titles.
        %         text(  ones( length(ctick),1 )*1.05,ctick+add,ctickstr, 'fontsize',fsize);
        %
        %
        %
        %     elseif nocbar==0 & icont==1 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
        %         clear ctickstr
        %         ctick=get(hc,'ytick');
        %         for j=1:length(ctick)
        %             nu=ctick(j);
        %             %te=strcat('10^','{',num2str(ctick(j)),'}' );
        %             if nu>99
        %                 %te=num2str(ctick(j),'%2.2e');
        %                 te=num2str(sigfig(ctick(j),2));
        %
        %             else
        %                % te=num2str(ctick(j),'%2.2e');
        %                 te=num2str(sigfig(ctick(j),2));
        %             end
        %             ctickstr(j,1:length(te))=te;
        %         end
        %      end


        %end

        if nocbar==0 & subplotting==1 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
            set(hc,'fontsize',fsize);
        elseif nocbar==0 & (bigcbar==0 | (bigcbar==1 & isub==length(idirs)) )
            set(hc,'fontsize',fsize-2);
        end




    end



    ylims_temp=get(h(isub).h,'ylim');   %re-scale to hide extra column put in to get the colorbars the same in both plots

    axis(h(isub).h,[timesTH(i).t(1) timesTH(i).t(pend(i)) ylims_temp]);


    if dan_test==1 %%% dan test %%%
        cax=get(gca,'clim');
        if iminovr==1
            cax(1)=mincovOvr;
        end
        if imaxovr==1
            cax(2)=maxcovOvr;
        end
        caxis(cax);
        %        hc=colorbar;
        %        set(hc,'fontsize',fsize-3);

    end


end  %for i=1:nplots2d


if clab==1 & icont==1 %round up contour values so labels match those on colorbar

    for iclabel_count=1:icont_neg+1

        if iclabel_count==1
            ch=cbfA(i).cpos;
            cb_dat = cbfA(i).cpos;
            c_handle = cbfB(i).cpos;
        else
            ch=cbfA(i).cneg;
            cb_dat = cbfA(i).cneg;
            c_handle = cbfB(i).cneg;
        end
        jc=1;
        while jc<size(cb_dat,2)
            %ch(1,jc)=str2num(num2str(10^(cbfA(i).c(1,jc)),'%2.2e'));

            if dlogflag==0 & logflag==1
                if 10^cb_dat(1,jc)>99
                    ch(1,jc)=sigfig(10^(cb_dat(1,jc) ),3);
                else
                    ch(1,jc)=sigfig(10^(cb_dat(1,jc) ),2);
                end
            elseif dlogflag==1
                %                    if idlog(cbfA(i).c(1,jc),dlogmin)>99
                %                         ch(1,jc)=sigfig(idlog(cbfA(i).c(1,jc),dlogmin ),2);
                %                     else
                %                         ch(1,jc)=sigfig(idlog(cbfA(i).c(1,jc),dlogmin ),2);
                %                     end
                ch(1,jc)=sigfig(idlog(cb_dat(1,jc),dlogmin),3);
            else
                ch(1,jc)=sigfig(cb_dat(1,jc),sig);
            end

            jc=jc+cb_dat(2,jc)+1;
        end

        if length(ch)>0
            if manclab==0
                clabel(ch,c_handle,'labelspacing',144); %default 144
            else
                clabel(ch,c_handle,'manual');
            end
        end

    end  %for iclabel_count=1:icont_neg+1
    icont_neg=0;

end %clab==1 & icont==1


%if clines==0
%   shading flat; %shading stops black contour lines from appearing
%end  %NOTE - needs to be after the contour labelling as this puts the
%lines back!
%This method produced incorrect results when have small NaN areas
%filled white - using contourf(.....,'linestyle','none') is better



if vectorf==1
    spy=25;
    spz=15;

    %    spy=20;
    %    spz=10;

    sqy=size(GridDan(idir).Y1(xinds),1)/spy;
    sqz=round((izmax-izmin)/spz);

    zinds=[izmin:sqz:izmax+2*sqz];
    yinds=[xinds(1):sqy:xinds(end)+sqy];

    sf=max(max(TwoD.V))/max(max(TwoD.W));
    hold on;
    quiver(GridDan(idir).Y1(yinds)./1000,GridDan(idir).Z(zinds)./1000,TwoD.V(zinds,yinds),TwoD.W(zinds,yinds),'w');
end


if itimestamp==1
    %    text(timesTH(i).t(1)-(timesTH(i).t(end)-timesTH(i).t(1))*0.12,((zz(1).z(end)-zz(1).z(1))*0.18+zz(1).z(end))/1000,['Time = ' timlab ' UTC'],'fontsize',18);
    %    text(timesTH(i).t(1)-(timesTH(i).t(end)-timesTH(i).t(1))*0.12,((zz(1).z(end)-zz(1).z(1))*0.08+zz(1).z(end))/1000,['Time = ' timlab ' UTC'],'fontsize',18);

    if subplotting==1
        %       text(0,0,['Time = ' timlab ' UTC'],'units','centimeters','position',[-1.5 16],'fontsize',fsize);
        %     text(0,0,['Time = ' timlab ' UTC'],'units','centimeters','position',[-1.5 13.2],'fontsize',fsize);
        %       text(0,0,['Time = ' timlab],'units','centimeters','position',[-1.5 13.2],'fontsize',fsize);
        text(0,0,['Time = ' timlab ' UTC'],'units','centimeters','position',[-1.5 18],'fontsize',fsize);

    else
        %  text(0,0,['Time = ' timlab ' UTC'],'units','centimeters','position',[-2.5 11.5],'fontsize',fsize);
        text(0,0,['Time = ' timlab],'units','centimeters','position',[-2.5 13],'fontsize',fsize);

    end


    % text(timesTH(i).t(1)*1.2,23.0,['Time = ' f ' UTC'],'fontsize',18);
end




%
% tims=[9:2:23];
% ti=datenum(2004,2,24,tims,0,0);
% set(gca,'xtick',[ti]);
% datetick('x',15,'keepticks');

if lememm==1   %rescale if doing lem/emm comparison so as to get same axes
    set(gca,'xlim',timelims);
    set(gca,'ylim',zlims);
    if length(clims)==2
        %        set(gca,'clim',clims);
    end
end


if idirstamp==1
    %     if subplotting==1
    %         ylims=get(h(iplot).h,'ylim');
    %         xlims=get(h(iplot).h,'xlim');
    %     else
    %         ylims=get(h.h,'ylim');
    %         xlims=get(h.h,'xlim');
    %     end

    ylims=get(h(iplot).h,'ylim');
    xlims=get(h(iplot).h,'xlim');

    %    text(timesTH(i).t(1)-0.5,ylims(2)*1.002,[direcDan(idir).dir],'fontsize',12);
    if subplotting==1
        dist=(ylims(2)-ylims(1))/7.8;
    else
        dist=(ylims(2)-ylims(1))/10.1;
    end

    %dist2=(timesTH(i).t(pend(i))-timesTH(i).t(1))/20;
    dist2=(xlims(end)-xlims(1))/8;
    %      dist2=(xlims(end)-xlims(1))/20; %last one used

    %    text(timesTH(i).t(1)-0.5,ylims(1)*0.988,[direcDan(idir).dir],'fontsize',12);
    axes(h(end).h);

    if plotcase==65
        dirname=run_name_emm_select;
    else
        dirname=runName(idir).nam
    end

    if iabc==1
        abc={'(a) ','(b) ','(c) ','(d) '};
        dirstr=[abc{iplot} dirname];
    else
        dirstr=dirname;
    end

    %    text(xlims(1)-dist2,ylims(1)-dist,[dirstr],'fontsize',fsize-2);
    text(xlims(1)-dist2,ylims(1)-dist,[dirstr],'fontsize',fsize-2);
    %	text(xlims(1)-4*dist2,ylims(2)+dist,[dirstr],'fontsize',fsize+4);
end


if (i2d~=1 & i2d~=2 & i2d~=3)
    xx=get(h(isub).h,'xticklabels');
    xx=str2num(xx);
    xx=num2str(mod(xx,24));
    set(h(isub).h,'xticklabels',xx);
end

set(gcf,'paperpositionmode','auto');


if isave==1
    set(gcf,'paperpositionmode','auto');
    print(gcf,'-djpeg','-r350',exname);
    %print(gcf,'-dmeta',exname);
    %close(gcf);
end



clims_terr = get(gca,'clim'); %need to remember and reset the color lims as goes weird when do contours

if iadd_terrain==1
    add_terrain;
end

if icolmap==1
    colormap(cmap);
end

if iplot_latlon==1
    plot_latlon_lines2;
end

if iadd_overlay==1
    hold on
    plot(x_overlay,y_overlay,'k');
end


caxis(clims_terr); %resest colourscale as sometimes contour causes it to go wrong

if iadd_wind_quivers==1

    [ax_orig]=plot_wind_quiver_arrows(u_quiver,v_quiver,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare);
    %for some reason it only plots the reference arrow head if run twice!
    %    plot_wind_quiver_arrows(u_quiver,v_quiver,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1);
    axes(ax_orig);  %revert back to the original axis for zooming in etc
    set(gca,'ylim',[zz(1).z(1) zz(1).z(end)]); %need to rescale the vertical axis
end

if iylim==1
    set(gca,'ylim',ylims);
end
if ixlim==1
    set(gca,'xlim',xlims);
end

if idraw_streamlines==1
    draw_streamlines;
end

%pause(60);


