function [Radar,ThreeD]=AFTERIMPORTDIAG7(ThreeD,Grid,RNc)% to be run after Importdiag
N_DROP=RNc;

% Calculate pressure
[r,c,p]=size(ThreeD.P);
PRESS=permute(repmat(Grid.PREFN,[1 c r]),[3 2 1])+ThreeD.P;
Tempera=ThreeD.TH2.*((PRESS)./100000).^0.286;

% Calculate density
RHO=(PRESS./287)./Tempera;

% moment parameters
[SM.NR0,SM.LAMDA_R]=SingleMom(ThreeD.Q(:,:,:,3),RHO,'rain');
[DM.NS0,DM.LAMDA_S]=DoubleMom(ThreeD.Q(:,:,:,9),ThreeD.Q(:,:,:,4),RHO,'snow');
[DM.NG0,DM.LAMDA_G]=DoubleMom(ThreeD.Q(:,:,:,8),ThreeD.Q(:,:,:,5),RHO,'graupel');
[DM.NI0,DM.LAMDA_I]=DoubleMom(ThreeD.Q(:,:,:,7),ThreeD.Q(:,:,:,6),RHO,'ice');
% radar parameters
DROP=(9.*ThreeD.Q(:,:,:,2)./2./RHO./1000./N_DROP./pi).^(1./3);
Radar.Zd=N_DROP.*(DROP./1E-3).^6;
[Radar.Zr]=real(ZfactorSingle(SM.NR0,SM.LAMDA_R,RHO,'rain'));
[Radar.Zs]=real(Zfactor(DM.NS0,DM.LAMDA_S,RHO,'snow'));
[Radar.Zg]=real(Zfactor(DM.NG0,DM.LAMDA_G,RHO,'graupel'));
[Radar.Zi]=real(Zfactor(DM.NI0,DM.LAMDA_I,RHO,'ice'));

% total 
II=find(isnan(Radar.Zr));Radar.Zr(II)=0.000000000000000000001;
II=find(isnan(Radar.Zs));Radar.Zs(II)=0.000000000000000000001;
II=find(isnan(Radar.Zg));Radar.Zg(II)=0.000000000000000000001;
II=find(isnan(Radar.Zi));Radar.Zi(II)=0.000000000000000000001;
II=find(isnan(Radar.Zd));Radar.Zd(II)=0.000000000000000000001;
clear II;
Radar.Z=Radar.Zr+0.224.*Radar.Zs+0.224.*Radar.Zg+0.224.*Radar.Zi+Radar.Zd;
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
ThreeD.RHO=RHO;
ThreeD.Temp=Tempera;
