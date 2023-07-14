%example plot of tephigram
groundheight=1000;

%[tp,h]=plot_tephi(-40+273,50+273,153,5000);

nz=length(GridDan(idir).Z);
pdat=GridDan(idir).PREFN(1:end);
qdat=GridDan(idir).OLQBAR(1:nz,1);
qdat(1)=qdat(2);
tdat=tempLES(GridDan(idir));
%tdat=tdat(2:end);
f=1e6*28.97/18;
qsat=satvappress(tdat,'goff','liq',pdat,1)/f;

%plot_tephi_data(tdat,pdat,qdat,qsat);
plot_tephi_data2(tdat,pdat,qdat,qsat,GridDan(idir).Z(1:end)+groundheight);