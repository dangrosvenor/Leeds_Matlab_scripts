%Put text values into the boxes of a pcolor plot% % Xmid_text = 0.5*( timesTH(1).t(1:end-1) + timesTH(1).t(2:end) );% Ymid_text = 0.5*( height(1:end-1) + height(2:end) );hold onfsize_txt=24;fac=1;%Find top 5 magnitudes and underline themdat = qh(1:end-1,1:end-1);  %Use the data from the 2D histogram, but ignoring the last row and columns since they are dummy values[Y,I]=sort(abs(dat(:)),1,'descend'); %Sort from largest to smallest[IY,IX]=ind2sub(size(dat),I);Ntop=10;%Remove NaNs since they get sorted to the top of the listinan=find(isnan(Y)==1);IY(inan)=[];IX(inan)=[];dx=diff(Xbins);for ix=1:length(mid_Xbins)    for iy=1:length(mid_Ybins)        top=0;        for itop=1:min(Ntop,length(IX))            if ix==IX(itop) & iy==IY(itop) %If it is one of the top N abs numbers                %text(mid_Xbins(ix)-dx(ix)/4,mid_Ybins(iy),['\textsf{\underline{' num2str(fac*qh(iy,ix),'%.2g') '}}'],'Interpreter','latex','fontsize',fsize_txt);                val = fac*qh(iy,ix);                if abs(val)<10                    prec='%.2f'; %2 decimal places                else                    prec='%.4g'; %4 sig figs max                end                text(mid_Xbins(ix)-dx(ix)/4,mid_Ybins(iy),['\textrm{{' num2str(val,prec) '}}'],'Interpreter','latex','fontsize',fsize_txt);                                %texsf uses the san serif font. The default is roman (if                %you don't specify, but also can use textrm. But this looks                %weird since it seems to miss the bottom of the                %characters.                top=1;            end        end        if top==0            %text(mid_Xbins(ix),mid_Ybins(iy),['\textsf{' num2str(fac*qh(iy,ix),'%.2g') '}'],'Interpreter','latex','fontsize',fsize_txt);                    end                    endend