% to be run after Importdiag
[YI,ZI]=meshgrid(Grid.Y1,Grid.Z);
[TwoD.PRESS,TwoD.RHO] = Press(TwoD.TH1',Grid.THREF',97570,YI,ZI);
% moment parameters
TwoD.PRESS=TwoD.PRESS';
TwoD.RHO=TwoD.RHO';
[SM.NR0,SM.LAMDA_R]=SingleMom(TwoD.Q(:,:,3),TwoD.RHO,'rain');
[DM.NS0,DM.LAMDA_S]=DoubleMom(TwoD.Q(:,:,9),TwoD.Q(:,:,4),TwoD.RHO,'snow');
[DM.NG0,DM.LAMDA_G]=DoubleMom(TwoD.Q(:,:,8),TwoD.Q(:,:,5),TwoD.RHO,'graupel');
[DM.NI0,DM.LAMDA_I]=DoubleMom(TwoD.Q(:,:,7),TwoD.Q(:,:,6),TwoD.RHO,'ice');
% radar parameters
[Radar.Zr]=real(ZfactorSingle(SM.NR0,SM.LAMDA_R,TwoD.RHO,'rain'));
[Radar.Zs]=real(Zfactor(DM.NS0,DM.LAMDA_S,TwoD.RHO,'snow'));
[Radar.Zg]=real(Zfactor(DM.NG0,DM.LAMDA_G,TwoD.RHO,'graupel'));
[Radar.Zi]=real(Zfactor(DM.NI0,DM.LAMDA_I,TwoD.RHO,'ice'));

% total 
II=find(isnan(Radar.Zr));Radar.Zr(II)=0.000000000000000000001;
II=find(isnan(Radar.Zs));Radar.Zs(II)=0.000000000000000000001;
II=find(isnan(Radar.Zg));Radar.Zg(II)=0.000000000000000000001;
II=find(isnan(Radar.Zi));Radar.Zi(II)=0.000000000000000000001;
clear II;
Radar.Z=Radar.Zr+Radar.Zs+Radar.Zg+Radar.Zi;

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
