	ip=findheight(GridDan(1).PREFN,450e2);
    
    TwoD=TwoDDan(3);
    
for k=1:size(TwoD.Q,1)
    iacc=find( sum(TwoD.Q(k,:,[10]),3) > 0.1 & TwoD.W(k,:)>0.5 & sum(TwoD.Q(k,:,[2 6]),3)>=1e-7  );    
 %   iacc=find(TwoD.Q(k,:,7)>1e6 );    
 %   iacc=find(sum(TwoD.Q(k,:,[2 6]),3)>1e-7 );    
 % iacc=find(sum(TwoD.Q(k,:,[6]),3)>1e-7 );    
 % iacc=find(TwoD.W(k,:)>1);  % & sum(TwoD.Q(k,:,[2 6]),3)>=1e-7);    

    P=TwoD.PP(k,iacc);
    pacc(k,:)=mean(P);
    T=TwoD.TH2(k,iacc)./(1e5./P).^0.286;
    Tacc(k,:)=mean(T);
    rho=P.*28.97e-3/8.3144./T;
    rhoacc(k)=mean(rho,2);
    LWCacc(k)=mean(sum(TwoD.Q(k,iacc,2),3).*rho);
    CONDacc(k)=mean(sum(TwoD.Q(k,iacc,2:6),3).*rho);
    RAINacc(k)=mean(sum(TwoD.Q(k,iacc,2:3),3).*rho);
    TRACERacc(k)=mean(sum(TwoD.Q(k,iacc,10),3).*rho);
    
    INCacc(k)=mean(sum(TwoD.Q(k,iacc,7:9),3).*rho);
    IWCacc(k)=mean(sum(TwoD.Q(k,iacc,4:6),3).*rho);
    
    ICEacc(k)=mean(sum(TwoD.Q(k,iacc,6),3).*rho);
    
    wind_acc(k)=mean(TwoD.W(k,iacc));


end

for k=1:size(TwoD.Q,1)
	iacc=find( sum(TwoD.Q(k,:,[4:6]),3) > 1e-4 );
    P=TwoD.PP(k,iacc);
    pacc2(k,:)=mean(P);
    T=TwoD.TH2(k,iacc)./(1e5./P).^0.286;
    Tacc2(k,:)=mean(T);
    rho=P.*28.97e-3/8.3144./T;
    rhoacc2(k)=mean(rho,2);
    LWCacc2(k)=mean(sum(TwoD.Q(k,iacc,2),3).*rho);
    CONDacc2(k)=mean(sum(TwoD.Q(k,iacc,2:6),3).*rho);
    RAINacc2(k)=mean(sum(TwoD.Q(k,iacc,2:3),3).*rho);
    TRACERacc2(k)=mean(sum(TwoD.Q(k,iacc,10),3).*rho);


end
'done EMM_ACCpressure.m'