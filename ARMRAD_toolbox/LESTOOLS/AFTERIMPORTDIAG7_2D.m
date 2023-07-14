function [Radar,TwoD]=AFTERIMPORTDIAG7_2D(TwoD,Grid,POINT,queryData,run2)% to be run after Importdiag
Q=find([queryData{:,7}]==run2);
N_DROP=queryData{Q,6};

% Calculate pressure
[r,c,p]=size(TwoD.P);
PRESS=permute(repmat(Grid.PREFN,[1 1 c]),[1 3 2])+TwoD.P;
Tempera=TwoD.TH2.*((PRESS)./100000).^0.286;

% Calculate density
RHO=(PRESS./287)./Tempera;

% moment parameters
[SM.NR0,SM.LAMDA_R]=SingleMom(squeeze(TwoD.Q(:,:,3,:)),RHO,'rain',queryData(Q,:));
[DM.NS0,DM.LAMDA_S]=DoubleMom(squeeze(TwoD.Q(:,:,9,:)),squeeze(TwoD.Q(:,:,4,:)),RHO,'snow',queryData(Q,:));
[DM.NG0,DM.LAMDA_G]=DoubleMom(squeeze(TwoD.Q(:,:,8,:)),squeeze(TwoD.Q(:,:,5,:)),RHO,'graupel',queryData(Q,:));
[DM.NI0,DM.LAMDA_I]=DoubleMom(squeeze(TwoD.Q(:,:,7,:)),squeeze(TwoD.Q(:,:,6,:)),RHO,'ice',queryData(Q,:));
% radar parameters
DROP=(9.*squeeze(TwoD.Q(:,:,2,:))./2./RHO./1000./N_DROP./pi).^(1./3);
Radar.Zd=N_DROP.*(DROP./1E-3).^6;
[Radar.Zr]=real(ZfactorSingle(SM.NR0,SM.LAMDA_R,RHO,'rain',queryData(Q,:)));
[Radar.Zs]=0.19.*real(Zfactor(DM.NS0,DM.LAMDA_S,RHO,'snow',queryData(Q,:)));
[Radar.Zg]=0.19.*real(ZfactorTripleMom(DM.NG0,DM.LAMDA_G,RHO,'graupel',POINT,TwoD,queryData,run2));% A CHANGER !!! faire une nouvelle fonction
[Radar.Zi]=0.19.*real(Zfactor(DM.NI0,DM.LAMDA_I,RHO,'ice',queryData(Q,:)));

% total 
II=find(isnan(Radar.Zr));Radar.Zr(II)=0.000000000000000000001;SM.NR0(II) = NaN;SM.LAMDA_R(II) = NaN;
II=find(isnan(Radar.Zs));Radar.Zs(II)=0.000000000000000000001;DM.NS0(II) = NaN;DM.LAMDA_S(II) = NaN;
II=find(isnan(Radar.Zg));Radar.Zg(II)=0.000000000000000000001;DM.NG0(II) = NaN;DM.LAMDA_G(II) = NaN;
II=find(isnan(Radar.Zi));Radar.Zi(II)=0.000000000000000000001;DM.NI0(II) = NaN;DM.LAMDA_I(II) = NaN;
II=find(isnan(Radar.Zd));Radar.Zd(II)=0.000000000000000000001;
clear II;
%on change le facteur qui apparait dans la ligne suivante, il etait 0.224
Radar.Z=Radar.Zr+Radar.Zs+Radar.Zg+Radar.Zi+Radar.Zd;
% Radar.Z=Radar.Zr+Radar.Zs+Radar.Zg+Radar.Zi+Radar.Zd;
%Radar.Z=Radar.Zr+0.224.*Radar.Zs+0.224.*Radar.Zi+Radar.Zd;
% set all nans to
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
%TwoD(M).Temp=(TwoD(M).TH1+TwoD(M).TH2).*(TwoD(M).PRESS./100000).^0.286;
TwoD.RHO=RHO;
TwoD.Press=PRESS;
TwoD.Temp=Tempera;
TwoD.SM.NR0=SM.NR0;
TwoD.SM.LAMDA_R=SM.LAMDA_R;
TwoD.DM.NS0=DM.NS0;
TwoD.DM.LAMDA_S=DM.LAMDA_S;
TwoD.DM.NG0=DM.NG0;
TwoD.DM.LAMDA_G=DM.LAMDA_G;
TwoD.DM.NI0=DM.NI0;
TwoD.DM.LAMDA_I=DM.LAMDA_I;