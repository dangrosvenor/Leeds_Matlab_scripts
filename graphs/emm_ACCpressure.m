ip=findheight(GridDan(1).PREFN,450e2);
for k=1:size(TwoD.Q,1)
    if k<ip
        iacc=find( sum(TwoD.Q(k,:,[2]),3) > 5e-4 );
    else
        iacc=find( sum(TwoD.Q(k,:,[2 4:6]),3) > 1e-4 );
    end
    
    P=TwoD.PP(k,iacc);
    pacc(k,:)=mean(P);
    T=TwoD.TH2(k,iacc)./(1e5./P).^0.286;
    Tacc(k,:)=mean(T);
    rho=P.*28.97e-3/8.3144./T;
    rhoacc(k)=mean(rho,2);
    LWCacc(k)=mean(sum(TwoD.Q(k,iacc,2),3).*rho);
    CONDacc(k)=mean(sum(TwoD.Q(k,iacc,2:6),3).*rho);
    RAINacc(k)=mean(sum(TwoD.Q(k,iacc,2:3),3).*rho);

end

'done EMM_ACCpressure.m'