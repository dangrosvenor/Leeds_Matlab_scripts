% to be run after Importdiag
[YI,ZI]=meshgrid(Grid(i).Y1,Grid(i).Z);
[TwoD(i).PRESS,TwoD(i).RHO] = Press(TwoD(i).TH,Grid(i).THREF,97670,YI,ZI);
% moment parameters
[SM(i).NR0,SM(i).LAMDA_R]=SingleMom(TwoD(i).Q(:,:,2),TwoD(i).RHO,'rain');
[DM(i).NS0,DM(i).LAMDA_S]=DoubleMom(TwoD(i).Q(:,:,8),TwoD(i).Q(:,:,3),TwoD(i).RHO,'snow');
[DM(i).NG0,DM(i).LAMDA_G]=DoubleMom(TwoD(i).Q(:,:,7),TwoD(i).Q(:,:,4),TwoD(i).RHO,'graupel');
[DM(i).NI0,DM(i).LAMDA_I]=DoubleMom(TwoD(i).Q(:,:,6),TwoD(i).Q(:,:,5),TwoD(i).RHO,'ice');
% radar parameters
[Radar(i).Zr]=ZfactorSingle(SM(i).NR0,SM(i).LAMDA_R,TwoD(i).RHO,'rain');
[Radar(i).Zs]=Zfactor(DM(i).NS0,DM(i).LAMDA_S,TwoD(i).RHO,'snow');
[Radar(i).Zg]=Zfactor(DM(i).NG0,DM(i).LAMDA_G,TwoD(i).RHO,'graupel');
[Radar(i).Zi]=Zfactor(DM(i).NI0,DM(i).LAMDA_I,TwoD(i).RHO,'ice');

% total 
Radar(i).Z=Radar(i).Zr+Radar(i).Zs+Radar(i).Zg+Radar(i).Zi;


find(Radar(i).Z < 10.^(-4));
Radar(i).Z(ans)=nan;
Radar(i).dbZ=10.*log10(Radar(i).Z);

% find temperature
TwoD(i).TEMP=(TwoD(i).PRESS./TwoD(i).RHO).*(28.*1.67E-27)./(1.381E-23);
TwoD(i).TEMP=TwoD(i).TEMP-273;
