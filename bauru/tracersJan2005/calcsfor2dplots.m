hrange=4;
switch hrange
case 1
    minZ=14e3;
    maxZ=23e3;
case 2
    minZ=10e3;
    maxZ=23e3;
case 3
    minZ=0.2e3;
    maxZ=23e3;
case 4
    minZ=14e3;
    maxZ=19e3; 
case 5
    minZ=13e3;
    maxZ=22e3;
end
%ncont=15;



z=GridDan(1).Z;

[izmin izmax]=findheight(z,minZ,maxZ);

switch i577
case 'temp_change_2d'
    P=GridDan(1).PREFN;      
    init=icediagsALL(1).i(:,1,246) ./(1e5./P).^0.286 ;  
    
    P=repmat(P,[1 length(Grid.Y1(2:end-1))]);
    T=TwoD.TH2(:,2:end-1)./(1000e2./P).^0.286;
    
    for iz=1:length(GridDan(idir).Z)
        [mx mi]=maxALL(TwoD.Q(iz,2:end-1,1));
        med=median( TwoD.Q(iz,2:end-1,1) );
        
        Tiz=T(iz,:);
        ivap=find( sum( TwoD.Q(iz,2:end-1,1:6) , 3) > med*1.05 ); %sum for total water
        dT_conv(1).dat(iz,jj) = squeeze( mean( Tiz(ivap)  ))  - init(iz);
        
        tracer_concs=[1e-4 1e-3 1e-2 1e-1];
        for itracer=1:length(tracer_concs)
            ivap=find( TwoD.Q(iz,2:end-1,10) >= tracer_concs(itracer) ); %sum for total water
            dT_conv_tracer(1).dat(iz,jj,itracer) = squeeze( mean( Tiz(ivap)  ))  - init(iz);
        end
        
        
        
        ivap=find( sum( TwoD.Q(iz,2:end-1,1:6) , 3) < med*1.05 ); %sum for total water
        dT_nonconv(1).dat(iz,jj) = squeeze( mean( Tiz(ivap)  ))  - init(iz);
        
        medT=median(median(Tiz));
        size_Tiz=size(Tiz);
        dT_bubble(1).dat(iz,jj) = squeeze( mean( Tiz([1:7 size_Tiz(2)-6:end]) ) )  - medT;
        
        
        tpertTimH_full(1).t(iz,:,:,jj)=Tiz( [size_Tiz(2)-6:end 1:7 ] ) - medT;
        Tpert_full(1).med=medT;
        
        tpert_max(1).t(iz,jj)=maxALL( Tiz( [size_Tiz(2)-6:end 1:7 ] ) - medT );
        
        vapiz=TwoD.Q(iz,2:end-1,1);
        medV=median(median(vapiz));          
        vappert_max(1).dat(iz,jj)=maxALL( vapiz( [size_Tiz(2)-6:end 1:7 ] ) - medV );
        vappert_full(1).dat(iz,:,:,jj)=vapiz( [size_Tiz(2)-6:end 1:7 ] ) - medV;
        vappert_full(1).med=medV;
        
    end

case 'temp_change_2d_2'
%    P=GridDan(1).PREFN;          
%    P=repmat(P,[1 length(Grid.Y1(2:end-1))]);
    P=TwoD.PP(:,2:end-1);
    T=TwoD.TH2(:,2:end-1)./(1000e2./P).^0.286;
    
    for iz=1:length(GridDan(idir).Z)
        %          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz) + totice_44(idir).dat(2:end-1,2:end-1,iz));
%         [mx mi]=maxALL(ThreeD.Q(2:end-1,2:end-1,iz));
%         med=median(median( ThreeD.Q(2:end-1,2:end-1,iz,1) ));
%         
         Tiz=T(iz,:);
%         ivap=find( sum( ThreeD.Q(2:end-1,2:end-1,iz,:) , 4) > med*1.05 ); %sum for total water
%         dT_conv(1).dat(iz,jj) = squeeze( mean( Tiz(ivap)  ))  - init(iz);
%         
%         tracer_concs=[1e-4 1e-3 1e-2 1e-1];
%         for itracer=1:length(tracer_concs)
%             ivap=find( ThreeD.Q(2:end-1,2:end-1,iz,2) >= tracer_concs(itracer) ); %sum for total water
%             dT_conv_tracer(1).dat(iz,jj,itracer) = squeeze( mean( Tiz(ivap)  ))  - init(iz);
%         end
        
        
%         
%         ivap=find( sum( ThreeD.Q(2:end-1,2:end-1,iz,:) , 4) < med*1.05 ); %sum for total water
%         dT_nonconv(1).dat(iz,jj) = squeeze( mean( Tiz(ivap)  ))  - init(iz);
        
        medT=median(median(Tiz));
         size_Tiz=size(Tiz);
%         dT_bubble(1).dat(iz,jj) = squeeze( mean(mean( Tiz(1,[1:7 size_Tiz(2)-6:end],[1:7 size_Tiz(3)-6:end]))  ))  - medT;
        
        
        tpertTimH_full(1).t(iz,:,jj)=Tiz(1, [ size_Tiz(2)-6:end 1:7] ) - medT;
        tpertTimH_full(1).med(iz,jj)=medT;
                
        tpert_max(1).t(iz,jj)=maxALL( Tiz(1, [ size_Tiz(2)-6:end 1:7] ) - medT );
        
        vapiz=TwoD.Q(iz,2:end-1,1);
        medV=median(median(vapiz));          
        vappert_max(1).dat(iz,jj)=maxALL( vapiz([ size_Tiz(2)-6:end 1:7] ) - medV );
        
        vappertTimH_full(1).t(iz,:,jj)=vapiz([ size_Tiz(2)-6:end 1:7] ) - medV;
        vappertTimH_full(1).med(iz,jj)=medV;
        
    end
    
case 'temp_change_3d'
    P=GridDan(1).PREFN;
%    init=icediagsALL(1).i(:,1,246) ./(1e5./P).^0.286 ;  
    [T]=temp_from_press_and_th(GridDan(idir),ThreeD.TH1,ThreeD.P);
    for iz=1:length(GridDan(idir).Z)
        %          [mx mi]=maxALL(ThreeDDan(idir).Q(2:end-1,2:end-1,iz) + totice_44(idir).dat(2:end-1,2:end-1,iz));
%         [mx mi]=maxALL(ThreeD.Q(2:end-1,2:end-1,iz));
%         med=median(median( ThreeD.Q(2:end-1,2:end-1,iz,1) ));
%         
         Tiz=T(iz,:,:);
%         ivap=find( sum( ThreeD.Q(2:end-1,2:end-1,iz,:) , 4) > med*1.05 ); %sum for total water
%         dT_conv(1).dat(iz,jj) = squeeze( mean( Tiz(ivap)  ))  - init(iz);
%         
%         tracer_concs=[1e-4 1e-3 1e-2 1e-1];
%         for itracer=1:length(tracer_concs)
%             ivap=find( ThreeD.Q(2:end-1,2:end-1,iz,2) >= tracer_concs(itracer) ); %sum for total water
%             dT_conv_tracer(1).dat(iz,jj,itracer) = squeeze( mean( Tiz(ivap)  ))  - init(iz);
%         end
        
        
%         
%         ivap=find( sum( ThreeD.Q(2:end-1,2:end-1,iz,:) , 4) < med*1.05 ); %sum for total water
%         dT_nonconv(1).dat(iz,jj) = squeeze( mean( Tiz(ivap)  ))  - init(iz);
        
        medT=median(median(Tiz));
         size_Tiz=size(Tiz);
%         dT_bubble(1).dat(iz,jj) = squeeze( mean(mean( Tiz(1,[1:7 size_Tiz(2)-6:end],[1:7 size_Tiz(3)-6:end]))  ))  - medT;
        
        
        tpertTimH_full(1).t(iz,:,jj)=Tiz(1, [ size_Tiz(2)-6:end 1:7] ) - medT;
        tpertTimH_full(1).med(iz,jj)=medT;
                
        tpert_max(1).t(iz,jj)=maxALL( Tiz(1, [ size_Tiz(2)-6:end 1:7],[size_Tiz(3)-6:end 1:7 ] ) - medT );
        
        vapiz=ThreeD.Q(2:end-1,2:end-1,iz,1);
        medV=median(median(vapiz));          
        vappert_max(1).dat(iz,jj)=maxALL( vapiz([ size_Tiz(2)-6:end 1:7],[size_Tiz(3)-6:end 1:7 ] ) - medV );
        
        vappertTimH_full(1).t(iz,:,:,jj)=vapiz([ size_Tiz(2)-6:end 1:7],[size_Tiz(3)-6:end 1:7 ] ) - medV;
        vappertTimH_full(1).med(iz,jj)=medV;
        
    end
    
    
case 'lwc_width'
    
    lwcs=[1e-7];  %threshold(s)
    
    for ilwcs=1:length(lwcs)
        [a b]=find(TwoD.Q(:,:,2)>=lwcs(ilwcs)); %all points with >=n dBz
        ihs=unique(a); %all the different height indices with points > threshold
        
        lwc_width(1).dat(1:length(Grid.Z),ilwcs,jj)=0;
        for iradar=1:length(ihs)
            lwc_width(1).dat(ihs(iradar),ilwcs,jj)=length(find(a==ihs(iradar))); %number of points > threshold at each height index
        end
    end
    
    
    
    
case 'supersat_vars'
    ih=findheight((Grid.Z+620)/1000 , 17.8);
    potemp_hslice(1).dat(:,:,jj)=ThreeD.TH1(:,:,ih);
    vap_hslice(1).dat(:,:,jj)=ThreeD.Q(:,:,ih);
    pressure_hslice(1).dat(:,:,jj)=ThreeD.P(:,:,ih);
    
    % snowMR_hslice(1).dat(:,:,jj)=ThreeD.Q(:,:,ih,1);
    %  iceMR_hslice(1).dat(:,:,jj)=ThreeD.Q(:,:,ih,1);    
    %  iceNC_hslice(1).dat(:,:,jj)=ThreeD.Q(:,:,ih,2);    
    % snowNC_hslice(1).dat(:,:,jj)=ThreeD.Q(:,:,ih,4);
    
case 'dqdehyd'
    'dq dehyd'
    f=1e6*28.97/18;
    
    totw=sum(TwoD.Q(:,:,[1:6]),3); %jc is the number of the eg. GridDan(jc) where the data is stored
    vap=sum(TwoD.Q(:,:,[1]),3);   %j just loops from 1:n
    
    j
    prcs=[0:5:100];
    tot_prctiles(1).t(1:size(totw,1),jj,1:length(prcs))=(prctile(totw',prcs))';
    vap_prctiles(1).t(1:size(vap,1),jj,1:length(prcs))=(prctile(vap',prcs))';
    
    minpps=[3.67 5 1 2 3 4]/f;
    %minpps=3.67/f;
    
    P=TwoD.PP; %tot P
    thref=repmat(Grid.THREF,[1 length(Grid.Y1)]);
    T=TwoD(idir).TH1+thref; %tot potemp
    T=T./(1e5./P).^0.286; %tot temp                                                
    rho=P.*28.97e-3/8.3144./T;
    rhomoist=rho.*(1+TwoD.Q(:,:,1))./(1+1.608*TwoD.Q(:,:,1));  
    
    for ipps=1:length(minpps)
        
        for ikm=1:size(totw,1)
            
            
            inon=find(totw(ikm,:)<minpps(ipps));
            dq_tot(1).d(ikm,jj,ipps)=sum(minpps(ipps)-totw(ikm,inon)) * f / size(totw,2);
            nn(1).n(ikm,jj,ipps)=length(inon);
            
            
            inon2=find(vap(ikm,:)<minpps(ipps));
            dq_vaps(1).d(ikm,jj,ipps)=sum(minpps(ipps)-vap(ikm,inon2)) * f / size(vap,2);
            nn2(1).n(ikm,jj,ipps)=length(inon2);
            
            
            
            dw = TwoD.W(ikm,:) - mean(TwoD.W(ikm,:));
            
            dvap=vap(ikm,:) - Grid.OLQBAR(ikm,1); %vapour pert
            
            iposup=find(dw > 0 & dvap > 0);
            inegup=find(dw > 0 & dvap < 0);
            
            iposd=find(dw < 0 & dvap > 0);
            inegd=find(dw < 0 & dvap < 0);
            
            vappos_up(1).dat(ikm,jj)=sum( dw(iposup) .* dvap(iposup) );
            vapneg_up(1).dat(ikm,jj)=sum( dw(inegup) .* dvap(inegup) );
            vappos_d(1).dat(ikm,jj)=sum( dw(iposd) .* dvap(iposd) );
            vapneg_d(1).dat(ikm,jj)=sum( dw(inegd) .* dvap(inegd) );
            
            dtot= totw(ikm,:) - sum(Grid.OLQBAR(ikm,[1:6]),2);
            
            iposup=find(dw > 0 & dtot > 0);
            inegup=find(dw > 0 & dtot < 0);
            
            iposd=find(dw < 0 & dtot > 0);
            inegd=find(dw < 0 & dtot < 0);
            
            totpos_up(1).dat(ikm,jj)=sum( dw(iposup) .* dtot(iposup) );
            totneg_up(1).dat(ikm,jj)=sum( dw(inegup) .* dtot(inegup) );
            totpos_d(1).dat(ikm,jj)=sum( dw(iposd) .* dtot(iposd) );
            totneg_d(1).dat(ikm,jj)=sum( dw(inegd) .* dtot(inegd) );
            
            
            
            if ipps==2
                rho_prof(1).tot(ikm,jj)=mean(rhomoist(ikm,inon),2); %mean density in low tot water points
                rho_prof(1).vap(ikm,jj)=mean(rhomoist(ikm,inon2),2); %same for low vapour points
                rho_prof(1).mean(ikm,jj)=mean(rhomoist(ikm,:),2); %mean for all points
            end
            
        end    
        
    end
    
    icediags_5thSept_2005_32;
    
    ih1=findheight(GridDan(idir).Z/1000+0.62,12);
    ih2=findheight(GridDan(idir).Z/1000+0.62,19);
    
    i3d=0;
    if i3d==0
        bins=[0:0.1:200];
        dat=sum(TwoD.Q(ih1:ih2,:,1:6),3)*f;
        totdist(1).v(:,jj)=binner(dat,bins);
        
        bins=[0:0.1:20];
        dat=TwoD.Q(ih1:ih2,:,1)*f;
        vapdist(1).v(:,jj)=binner(dat,bins);
    else               
        bins=[0:0.1:30];
        %dat=ThreeD.Q(2:end-1,:,ih1:ih2)*f;
        for ik=1:ih2-ih1+1
            vapdist(1).v(:,jj,ik)=binner(f*ThreeD(1).Q(2:end-1,:,ik+ih1-1),bins); 
            meanvap(1).m(:,jj,ik)=squeeze(mean(mean(ThreeD(1).Q(2:end-1,:,ik+ih1-1)) ));
            
        end
    end
    
case 'cape'
    cape_vals;
    
case 'w_3d'
    MaxW.w(:,jj)=squeeze(max(max(ThreeD.W(:,:,:),[],2),[],1)); %profiles of max w
    MinW.w(:,jj)=squeeze(min(min(ThreeD.W(:,:,:),[],2),[],1)); %profiles of max w
    
case 'vertvel'
    getLEM_up_width;
    
case 'meaninc'
    for k=1:size(TwoD.Q,1)
        iacc=find( sum(TwoD.Q(k,:,[7]),3) > 1e8 );  %& TwoD.W(k,:)>2 );    
        %  iacc=find(TwoD.W(k,:)>5 );    
        
        P2=TwoD.PP(k,:);
        P=TwoD.PP(k,iacc);
        pacc(k,:)=mean(P);
        T=TwoD.TH2(k,iacc)./(1e5./P).^0.286;
        T2=TwoD.TH2(k,:)./(1e5./P2).^0.286;
        Tacc(k,:)=mean(T);
        rho=P.*28.97e-3/8.3144./T;
        rho2=P2.*28.97e-3/8.3144./T2;
        rhoacc_alltim(1).dat(k,jj)=mean(rho,2);
        INCacc_alltim(1).dat(k,jj)=mean(sum(TwoD.Q(k,iacc,7:9),3).*rho);
        INCmaxacc_alltim(1).dat(k,jj)=max(sum(TwoD.Q(k,:,7:9),3).*rho2,[],2);    
        IWCacc_alltim(1).dat(k,jj)=mean(sum(TwoD.Q(k,iacc,4:6),3).*rho);
        
    end
    
    
case 'tracer'
    for k=1:size(TwoD.Q,1)
        iacc=find( sum(TwoD.Q(k,:,[10]),3) > 1e-7 & TwoD.W(k,:)>2 ); 
        iacc=find( sum(TwoD.Q(k,:,[10]),3) > 0.1 & TwoD.W(k,:)>0.5 & sum(TwoD.Q(k,:,[2 6]),3)>=1e-7  );    
        
        iacc=1:length(Grid.Y1);
        
        P=TwoD.PP(k,iacc);
        pacc(k,:)=mean(P);
        T=TwoD.TH2(k,iacc)./(1e5./P).^0.286;
        Tacc(k,:)=mean(T);
        rho=P.*28.97e-3/8.3144./T;
        rhoacc3(1).dat(k,jj)=mean(rho,2);
        LWCacc3(1).dat(k,jj)=mean(sum(TwoD.Q(k,iacc,2),3).*rho);
        CONDacc3(1).dat(k,jj)=mean(sum(TwoD.Q(k,iacc,2:6),3).*rho);
        RAINacc3(1).dat(k,jj)=mean(sum(TwoD.Q(k,iacc,2:3),3).*rho);
        TRACERacc3(1).dat(k,jj)=mean(sum(TwoD.Q(k,iacc,10),3).*rho);
        TRACERmax(1).dat(k,jj)=max(TwoD.Q(k,:,10).*rho,[],2);        
    end
    
case 'temppert_signchange_timser_LNB'
    TempPertChangeSign
    tpertsign_lnb(1).t(jj)=h; %estimate of LNB from where temp pert changes sign
    
case 'vapdists'
    ih1=findheight(GridDan(idir).Z/1000+0.62,12);
    ih2=findheight(GridDan(idir).Z/1000+0.62,19);
    
    i3d=0;
    if i3d==0
        bins=[0:0.1:200];
        for ik=1:ih2-ih1+1
            totdist(1).v(:,jj,ik)=binner(sum(f*TwoD.Q(ik+ih1-1,:,1:6),3),bins);
        end
        
        bins=[0:0.1:20];
        
        for ik=1:ih2-ih1+1
            vapdist(1).v(:,jj,ik)=binner(f*TwoD.Q(ik+ih1-1,:,1),bins);
        end
    else               
        bins=[0:0.1:30];
        %dat=ThreeD.Q(2:end-1,:,ih1:ih2)*f;
        for ik=1:ih2-ih1+1
            vapdist(1).v(:,jj,ik)=binner(f*ThreeD(1).Q(2:end-1,:,ik+ih1-1),bins); 
            meanvap(1).m(:,jj,ik)=squeeze(mean(mean(ThreeD(1).Q(2:end-1,:,ik+ih1-1)) ));
            
        end
    end
    
    
    
case 'winds'
    TwoD_alltim.W(:,:,jj)=TwoD(1).W;
    TwoD_alltim.V(:,:,jj)=TwoD(1).V;
    
    
case 'radar'
    %            ztot=Radar(Grid,TwoD,1,size(TwoD.Q,1));
    
    [r,c,p]=size(TwoD.P);
    PRESS=permute(repmat(Grid.PREFN,[1 1 c]),[1 3 2])+TwoD.P;
    Tempera=TwoD.TH2.*((PRESS)./100000).^0.286;
    % 
    % Calculate density
    RHO=(PRESS./287)./Tempera;
    
    ztot=Radar_new(GridDan(idir),TwoD,1,size(TwoD.Q,1),RHO);
    
    ztot=10*log10(ztot); 
    
    rads=[30 20 15 10 40 35];  %was [30 20 15 10] until 8th Jan, 2007
    
    for irad=1:length(rads)
        [a b]=find(ztot>=rads(irad)); %all points with >=n dBz
        ihs=unique(a); %all the different height indices with points > 10 dBz
        
        n10dbz(1).n(1:length(Grid.Z),irad,jj)=0;
        for iradar=1:length(ihs)
            n10dbz(1).n(ihs(iradar),irad,jj)=length(find(a==ihs(iradar))); %number of points > n dBz at each height index with > 10 dBz
        end
    end
    
    [maxOverall imax]=maxALL(ztot);
    n10dbz(1).max(jj)=maxOverall;
    n10dbz(1).hmax(jj)=Grid.Z(imax(1))+620;
    
case 'temps'
    
    if  jj==fnmin %only search for column numbers on first pass 
        clear dgcol
        id=1;
        dgcol(jfile).d(id)=getDGcol('ALL_TH',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ACC_TH',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALu_TH',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALd_TH',dgstrDan(jc).dg);
        id=id+1;
    end
    
    for idg=1:length(dgcol(jfile).d)
        icediagsTEMP(jfile).i(:,jj,idg)=TimeAvDan(jc).DGAV(:,dgcol(jfile).d(idg)); %saved as icediagsALL_a-b from Nov 2005
    end
    
case 'lnb'       
    LNB_2d; %do new LNB calculations            
    
    % 			pdat(1).p=lnb2d(izmin:izmax,:);  
    %             
    %             minlnb(jfile).m(:,jj)=min(lnb2d,[],2);
    %             maxlnb(jfile).m(:,jj)=max(lnb2d,[],2);
    %             
    %             zref=repmat(GridDan(1).Z(1:size(lnb2d,1))/1000+add_ground_height,[1 size(lnb2d,2)]);
    %             lnbdiff=lnb2d-zref;
    %             meanlnb_abv(jfile).m(:,jj)=zref(:,1)+meanselect(lnbdiff,'dat>0'); %calculate the mean only for points where lnb is lower than where air at
    %             meanlnb_bel(jfile).m(:,jj)=zref(:,1)+meanselect(lnbdiff,'dat<0'); %calculate the mean only for points where lnb is higher than where air at            
    %             
    %             bins(jfile).b=GridDan(1).Z/1000+add_ground_height;
    %             ipos=find(lnbdiff>=0);
    %             ineg=find(lnbdiff<0);
    %             
    %             lnbtemp=lnb2d;
    %             lnbtemp(ineg)=0;
    %             lnbbins_pos(jfile).l(:,jj)=binner(lnbtemp,bins(jfile).b); %put lnbs into bins - positiviely buoyant only
    %             
    %             lnbtemp=lnb2d;
    %             lnbtemp(ipos)=0;
    %             lnbbins_neg(jfile).l(:,jj)=binner(lnbtemp,bins(jfile).b); %put lnbs into bins - negatively buoyant only
    
    if  jj==fnmin 
        clear minlnb_vap maxlnb_vap meanlnb_abv_vap meanlnb_bel_vap lnbbins_pos_vap lnbbins_neg_vap bins_vap
        clear minlnb_tot maxlnb_tot meanlnb_abv_tot meanlnb_bel_tot lnbbins_pos_tot lnbbins_neg_tot bins_tot
    end
    
    [minlnb_vap(jfile).m(:,jj),maxlnb_vap(jfile).m(:,jj),meanlnb_abv_vap(jfile).m(:,jj),meanlnb_bel_vap(jfile).m(:,jj),lnbbins_pos_vap(jfile).l(:,jj),lnbbins_neg_vap(jfile).l(:,jj),bins_vap]...
        =lnb_calcs(lnb2d_vap,GridDan,add_ground_height);
    [minlnb_tot(jfile).m(:,jj),maxlnb_tot(jfile).m(:,jj),meanlnb_abv_tot(jfile).m(:,jj),meanlnb_bel_tot(jfile).m(:,jj),lnbbins_pos_tot(jfile).l(:,jj),lnbbins_neg_tot(jfile).l(:,jj),bins_tot]...
        =lnb_calcs(lnb2d_tot,GridDan,add_ground_height);
    
    if  jj==fnmin %only search for column numbers on first pass 
        clear dgcol
        id=1;
        dgcol(jfile).d(id)=getDGcol('ALL_LW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALL_SW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ACC_LW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ACC_SW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALu_LW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALu_SW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALd_LW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALd_SW',dgstrDan(jc).dg);
        id=id+1;
    end
    
    for idg=1:length(dgcol(jfile).d)
        icediagsRAD(jfile).i(:,jj,idg)=TimeAvDan(jc).DGAV(:,dgcol(jfile).d(idg)); %saved as icediagsALL_a-b from Nov 2005
    end
    
case 'rad'
    if  jj==fnmin %only search for column numbers on first pass 
        clear dgcol
        id=1;
        dgcol(jfile).d(id)=getDGcol('ALL_LW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALL_SW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ACC_LW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ACC_SW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALu_LW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALu_SW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALd_LW',dgstrDan(jc).dg);
        id=id+1;
        dgcol(jfile).d(id)=getDGcol('ALd_SW',dgstrDan(jc).dg);
        id=id+1;
    end
    
    for idg=1:length(dgcol(jfile).d)
        icediagsRAD(jfile).i(:,jj,idg)=TimeAvDan(jc).DGAV(:,dgcol(jfile).d(idg)); %saved as icediagsALL_a-b from Nov 2005
    end
    
    
case 'fall'
    ix=[1:length(GridDan(jc).Y1)];
    for iz=100:170 %length(GridDan(jc).Z);  
        rho=GridDan(jc).RHON(iz);
        for it=1:3
            switch it
            case 1 %ice
                lab='Ice';
                im=6;
                in=7;
            case 2 %snow
                lab='Snow';
                im=4;
                in=9;
            case 3 %graupel
                lab='Graupel';
                im=5;
                in=8;
            end
            
            q=TwoD.Q(iz,ix,im);
            n=TwoD.Q(iz,ix,in);
            [lam,nx0,dist(jfile).n(:,iz,jj,it),dist(jfile).m(:,iz,jj,it),D]=gamlemRow(n,q,rho,lower(lab(1)));
        end
    end
    
case 'ozone'
    pdat(1).p=TwoD(idir).Q(izmin:izmax,:,15); %vertical velocity    
case 'potemp'
    tref=repmat(GridDan(idir).THREF(izmin:izmax),[1 length(GridDan(idir).Y1)]);
    pdat(1).p=squeeze(sum(TwoD(idir).TH1(izmin:izmax,:),3))+tref; %potemp
case 'wind'
    pdat(1).p=squeeze(sum(TwoD(idir).W(izmin:izmax,:),3)); %vertical velocity
case 'lowtracer'                
    pdat(1).p=squeeze(sum(TwoD(idir).Q(izmin:izmax,:,[10]),3)); %low tracer
    maxtra(1).m(:,jj)=max(TwoD(idir).Q(:,:,[10]),[],2); %low tracer
    
case 'totwater'                
    pdat(1).p=f*squeeze(sum(TwoD(idir).Q(izmin:izmax,:,[1:6]),3)); %total water
case 'vapour'                
    pdat(1).p=f*squeeze(sum(TwoD(idir).Q(izmin:izmax,:,[1]),3)); %vapour
case 'tot_condensate'                
    %            pdat(1).p=1000*squeeze(sum(TwoD(idir).Q(izmin:izmax,:,[2:6]),3)); %tot condensate
    pdat(1).p=1000*squeeze(sum(TwoD.Q(izmin:izmax,:,[2:6]),3)); %tot condensate
case 'si'
    tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
    T=squeeze(sum(TwoD(idir).TH1,3))+tref;
    P=TwoD(idir).PP;
    T=T./(1e5./P).^0.286;
    qsi=satvapPress(T,'lem','ice',P,1)/fact; %satvappress gives in ppmv if 5th argument=1
    si=100*(TwoD(idir).Q(:,:,1)-qsi)./qsi;
    pdat(1).p=si(izmin:izmax,:);
    
    simaxTimH(1).s(:,jj)=max(si,[],2);
    siminTimH(1).s(:,jj)=min(si,[],2);
    simean(1).s(:,jj)=mean(si,2);
    
case 'temppert'
    tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]); %ref potemp
    T=squeeze(sum(TwoD(idir).TH1,3))+tref; %potemp
    pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
    P=TwoD(idir).PP; %actual p
    T=T./(1e5./P).^0.286; %actual T
    tref=tref./(1e5./pref).^0.286; %ref temp
    
    Tp=T-tref; %perturbation of temperature
    
    pdat(1).p=Tp(izmin:izmax,:);
    
    tpertTimH(1).t(:,jj)=mean(Tp,2); %mean temp pert
    
    tpertTimH_full(1).t(:,:,jj)=Tp;
    
case 'vapfull'
    'vapfull'
    vapfull(1).v(:,:,jj)=TwoD.Q(:,:,1);
    
    
    
case 'rhopert'
    tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
    pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
    
    
    th=TwoD(idir).TH1+tref; %tot potemp
    Tav=repmat(mean(th,2),[1 length(GridDan(idir).Y1)]); %mean T at this point in time
    
    P=TwoD(idir).PP; %tot P
    Pav=repmat(mean(P,2),[1 length(GridDan(idir).Y1)]); %mean P at this time
    
    T=th./(1e5./P).^0.286; %tot temp
    tref=tref./(1e5./pref).^0.286; %tot temp
    
    Tav=Tav./(1e5./Pav).^0.286; %tot temp
    
    rho=P.*28.97e-3/8.3144./T;
    
    
    %rhoref=pref.*28.97e-3/8.3144./tref;
    rhoav=Pav.*28.97e-3/8.3144./Tav; %average density at each time
    rhoav=rhoav.*(1+TwoD(idir).Q(:,:,1))./(1+1.608*TwoD(idir).Q(:,:,1)); %moist average
    
    if jj==1 %if is the inital run then make this the reference so are calculating changes from inital state *** need to run from run 1 ***
        rhoref=rhoav;
    end
    
    rhomoist=rho.*(1+TwoD(idir).Q(:,:,1))./(1+1.608*TwoD(idir).Q(:,:,1));
    
    rhopert=rhomoist-rhoref;
    rhopert2=rhomoist-rhoav;
    
    pdat(1).p=rhopert(izmin:izmax,:);
    
    rhopertTimH(1).t(:,jj)=mean(rhopert,2); %mean temp pert
    rhopertTimHmax(1).t(:,jj)=max(rhopert,[],2); %mean temp pert
    rhopertTimHmin(1).t(:,jj)=min(rhopert,[],2); %mean temp pert
    
    rhopertTimH2(1).t(:,jj)=mean(rhopert2,2); %mean temp pert
    rhopertTimHmax2(1).t(:,jj)=max(rhopert2,[],2); %mean temp pert
    rhopertTimHmin2(1).t(:,jj)=min(rhopert2,[],2); %mean temp pert
    
    
    totw=sum(TwoD.Q(:,:,[1:6]),3); %jc is the number of the eg. GridDan(jc) where the data is stored
    vap=sum(TwoD.Q(:,:,[1]),3);   %j just loops from 1:n
    prcs=[0:5:100];
    
    for ikm=100:200
        inon=find(vap(ikm,:)<5/f); %find indices for points with <5 ppmv
        inon2=find(totw(ikm,:)<5/f);
        rho5ppmv_vap(1).r(ikm,jj)=mean(rhopert(ikm,inon));
        rho5ppmv_tot(1).r(ikm,jj)=mean(rhopert(ikm,inon2));
        
        inon=find(vap(ikm,:)<5/f & rhopert(ikm,:)<0); %find indices for points with <5 ppmv and neg buoyant 
        inon2=find(vap(ikm,:)<5/f & rhopert(ikm,:)>=0); %pos buoyant
        rho5ppmv_vapneg(1).r(ikm,jj)=mean(rhopert(ikm,inon)); %mean density change
        rho5ppmv_vappos(1).r(ikm,jj)=mean(rhopert(ikm,inon2));
        
        inon=find(totw(ikm,:)<5/f & rhopert(ikm,:)<0); %find indices for points with <5 ppmv
        inon2=find(totw(ikm,:)<5/f & rhopert(ikm,:)>=0);
        rho5ppmv_totneg(1).r(ikm,jj)=mean(rhopert(ikm,inon));
        rho5ppmv_totpos(1).r(ikm,jj)=mean(rhopert(ikm,inon2));
        
        
        potemp_vap(1).r(ikm,jj)=mean(th(ikm,inon));
        potemp_tot(1).r(ikm,jj)=mean(th(ikm,inon2));
        
        for ipr=1:length(prcs)
            for itr=10:10
                iprf=find(TwoD.Q(ikm,:,itr)>prcs(ipr)/100); %index for all points with tracer val larger than prcs/100 (max = 1 g/kg)
                tra(1).mass(ikm,jj,itr,ipr)=sum(TwoD.Q(ikm,iprf,itr)); %total mass of tracer above different vals
                tra(1).num(ikm,jj,itr,ipr)=length(iprf); %store number of points as well                        
            end
        end
        
    end                
    
    
    
    low_prctiles(1).t(1:size(TwoD.Q(:,:,10),1),jj,1:length(prcs))=(prctile(TwoD.Q(:,:,10)',prcs))';
    % mid_prctiles(1).t(1:size(TwoD.Q(:,:,11),1),jj,1:length(prcs))=(prctile(TwoD.Q(:,:,11)',prcs))';
    % upp_prctiles(1).t(1:size(TwoD.Q(:,:,12),1),jj,1:length(prcs))=(prctile(TwoD.Q(:,:,12)',prcs))';
    % ztr_prctiles(1).t(1:size(TwoD.Q(:,:,13),1),jj,1:length(prcs))=(prctile(TwoD.Q(:,:,13)',prcs))';
    
    
    
case 'hydbal'
    %pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref presssure
    % pdat(1).p=TwoD(idir).PP(izmin:izmax,:)-pref(izmin:izmax,:);
    %zrefdiff=repmat(diff(GridDan(idir).Z(izmin-1:izmax)),[1 length(GridDan(idir).Y1)]); %height diff
    zref=repmat(GridDan(idir).Z(izmin-1:izmax+1),[1 length(GridDan(idir).Y1)]);
    
    %dpdz=diff(TwoD(idir).PP(izmin-1:izmax,:),1,1)./zrefdiff;
    
    dpdz=diffdan(TwoD(idir).PP(izmin-1:izmax+1,:),zref,1);
    
    tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
    T=TwoD(idir).TH1+tref; %tot potemp
    P=TwoD(idir).PP; %tot P
    T=T./(1e5./P).^0.286; %tot temp
    rho=P.*28.97e-3/8.3144./T;
    
    
    rhomoist=rho.*(1+TwoD(idir).Q(:,:,1))./(1+1.608*TwoD(idir).Q(:,:,1));
    
    pdat(1).p=-rhomoist(izmin:izmax,:)*9.81 - dpdz; %hydrostatic balance : remaining upwards force
    
case 'rhog'
    tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
    pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
    
    
    T=TwoD(idir).TH1+tref; %tot potemp
    Tav=repmat(mean(T,2),[1 length(GridDan(idir).Y1)]); %mean T at this point in time
    
    P=TwoD(idir).PP; %tot P
    Pav=repmat(mean(P,2),[1 length(GridDan(idir).Y1)]); %mean P at this time
    
    T=T./(1e5./P).^0.286; %tot temp
    tref=tref./(1e5./pref).^0.286; %tot temp
    
    Tav=Tav./(1e5./Pav).^0.286; %tot temp
    
    rho=P.*28.97e-3/8.3144./T;
    rhoref=pref.*28.97e-3/8.3144./tref;
    %rhoref=Pav.*28.97e-3/8.3144./Tav;
    
    
    rhomoist=rho.*(1+TwoD(idir).Q(:,:,1))./(1+1.608*TwoD(idir).Q(:,:,1));
    
    rhopert=rhomoist-rhoref;
    
    pdat(1).p=9.81*rhopert(izmin:izmax,:);
    
    %            rhopertTimH(1).t(:,jj)=mean(rhopert,2); %mean temp pert
    
    
case 'dpdz'
    %pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref presssure
    % pdat(1).p=TwoD(idir).PP(izmin:izmax,:)-pref(izmin:izmax,:);
    %zrefdiff=repmat(diff(GridDan(idir).Z(izmin-1:izmax)),[1 length(GridDan(idir).Y1)]); %height diff
    zref=repmat(GridDan(idir).Z(izmin-1:izmax+1),[1 length(GridDan(idir).Y1)]);
    
    %dpdz=diff(TwoD(idir).PP(izmin-1:izmax,:),1,1)./zrefdiff;
    
    dpdz=diffdan(TwoD(idir).PP(izmin-1:izmax+1,:),zref,1);
    
    tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
    T=TwoD(idir).TH1+tref; %tot potemp
    P=TwoD(idir).PP; %tot P
    T=T./(1e5./P).^0.286; %tot temp
    rho=P.*28.97e-3/8.3144./T;
    
    dpdz=diffdan(TwoD(idir).PP(izmin-1:izmax+1,:),zref,1);
    
    pref=repmat(GridDan(idir).PREFN,[1 length(GridDan(idir).Y1)]); %ref p
    refdpdz=diffdan(pref(izmin-1:izmax+1,:),zref,1);
    
    pdat(1).p= - dpdz + refdpdz; %pressure grad perturbation
    
    
    
case 'radiation'
    pdat(1).p=TwoD(idir).twdrad(izmin:izmax,:);
    
    
    
    
end    


