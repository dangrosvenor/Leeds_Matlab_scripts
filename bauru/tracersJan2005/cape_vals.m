idir=1;

%T=squeeze(sum(TwoD(idir).TH2,3)); %potemp
%P=TwoD(idir).PP; %actual p
%T=T./(1e5./P).^0.286; %actual T
%R=TwoD.Q(:,:,1);
ALT=GridDan(1).Z;

%qsat=satvappress(T,'goff','liq',P,1)/f;

D=7; %km
dx=0.25;
iz=2;
ix1=D/dx;
ix2=4098-D/dx;

% Tmean=mean([T(iz:end,1:ix1) T(iz:end,ix2:end)],2);
% Pmean=mean([P(iz:end,1:ix1) P(iz:end,ix2:end)],2);
% Rmean=mean([R(iz:end,1:ix1) R(iz:end,ix2:end)],2);
% QSATmean=mean([qsat(iz:end,1:ix1) qsat(iz:end,ix2:end)],2);

f=1e6*28.97/18;

ih=findheight(GridDan(1).Z,2.5e3);
Tback=temples(GridDan(1));
Pmean=GridDan(1).PREFN(2:end);

for ix=1:10
    Tmean=squeeze(mean7km(2:end,:,ix)) + Tback(2:end);
    Rmean=squeeze(vmean7km(2:end,:,ix))  + backgroundQ(2:end,:,ix) ;
    QSATmean=satvappress(Tmean,'goff','liq',Pmean,1)/f;

    
    
	[CAPE(ix),CIN(ix),HLCL(ix),TLCL(ix)]=CALC_CAPE(squeeze(Pmean),squeeze(Tmean),squeeze(Rmean),squeeze(QSATmean),ALT(2:end));
end

break


for ix=1:100  %size(TwoD.PP,2)
%    [CAPE(ix),CIN(ix),HLCL(ix),TLCL(ix)]=CALC_CAPE_v1(squeeze(P(iz:end,ix)),squeeze(T(iz:end,ix)),squeeze(R(iz:end,ix)),ALT(iz:end));
        [CAPE(ix),CIN(ix),HLCL(ix),TLCL(ix)]=CALC_CAPE(squeeze(P(iz:end,ix)),squeeze(T(iz:end,ix)),squeeze(R(iz:end,ix)),squeeze(qsat(iz:end,ix)),ALT(iz:end));

end

cape_etc(1).CAPE(:,jj)=CAPE;
cape_etc(1).CIN(:,jj)=CIN;
cape_etc(1).HLCL(:,jj)=HLCL;
cape_etc(1).TLCL(:,jj)=TLCL;

'done Cape calc'