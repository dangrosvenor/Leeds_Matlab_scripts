function [Radar,TwoD]=AFTERIMPORTDIAG6(TwoD,Grid,RNc,M)% to be run after Importdiag
N_DROP=RNc;
[YI,ZI]=meshgrid(Grid(M).Y1,Grid(M).Z);
[TwoD(M).PRESS,TwoD(M).RHO] = Press(TwoD(M).TH1',Grid(M).THREF',Grid(M).PREFN(1),YI,ZI);
% moment parameters
TwoD(M).PRESS=TwoD(M).PRESS';
TwoD(M).RHO=TwoD(M).RHO';
[SM.NR0,SM.LAMDA_R]=SingleMom(TwoD(M).Q(:,:,3),TwoD(M).RHO,'rain');
[DM.NS0,DM.LAMDA_S]=DoubleMom(TwoD(M).Q(:,:,9),TwoD(M).Q(:,:,4),TwoD(M).RHO,'snow');
[DM.NG0,DM.LAMDA_G]=DoubleMom(TwoD(M).Q(:,:,8),TwoD(M).Q(:,:,5),TwoD(M).RHO,'graupel');
[DM.NI0,DM.LAMDA_I]=DoubleMom(TwoD(M).Q(:,:,7),TwoD(M).Q(:,:,6),TwoD(M).RHO,'ice');
% radar parameters
DROP=(9.*TwoD(M).Q(:,:,2)./2./TwoD(M).RHO./1000./N_DROP./pi).^(1./3);
Radar.Zd=N_DROP.*(DROP./1E-3).^6;
[Radar.Zr]=real(ZfactorSingle(SM.NR0,SM.LAMDA_R,TwoD(M).RHO,'rain'));
[Radar.Zs]=real(Zfactor(DM.NS0,DM.LAMDA_S,TwoD(M).RHO,'snow'));
[Radar.Zg]=real(Zfactor(DM.NG0,DM.LAMDA_G,TwoD(M).RHO,'graupel'));
[Radar.Zi]=real(Zfactor(DM.NI0,DM.LAMDA_I,TwoD(M).RHO,'ice'));

% total 
II=find(isnan(Radar.Zr));Radar.Zr(II)=0.000000000000000000001;
II=find(isnan(Radar.Zs));Radar.Zs(II)=0.000000000000000000001;
II=find(isnan(Radar.Zg));Radar.Zg(II)=0.000000000000000000001;
II=find(isnan(Radar.Zi));Radar.Zi(II)=0.000000000000000000001;
II=find(isnan(Radar.Zd));Radar.Zd(II)=0.000000000000000000001;
clear II;
Radar.Z=Radar.Zr+0.224.*Radar.Zs+0.224.*Radar.Zg+0.224.*Radar.Zi+Radar.Zd;

size(Radar.Z);
jmax=ans(1);
size(Radar.Z);
kmax=ans(2);
%for j=1:jmax
%    for k=1:kmax
%        if((Radar.Z(j,k)) < 10.^(-4))
            %Zx(j,k)=10.^(-100./10.);
%            Radar.Z(j,k)=nan;
%end
%end
%end
TwoD(M).Temp=(TwoD(M).TH1+TwoD(M).TH2).*(TwoD(M).PRESS./100000).^0.286;