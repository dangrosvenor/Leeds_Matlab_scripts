hmtype='snownc';

switch hmtype
case 'ice'
dgs{1}='PIMLT';
dgs{2}='PSAUT';
dgs{3}='PSACI';
dgs{4}='PRACI_S';
dgs{5}='PGACI';
dgs{6}='PRACI_G';
dgs{7}='PIHAL';
dgs{8}='PIPRM';
dgs{9}='PICNT';
dgs{10}='PIDEP';
dgs{11}='PIACW';
dgs{12}='PIFRW';

case 'snownc'
riacr_s=getDGAVs('ALL_RIACR_S',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
rsaut=getDGAVs('ALL_RSAUT',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
rsbrk=getDGAVs('ALL_RSBRK',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
pssub=getDGAVs('ALL_PSSUB',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
psmlt=getDGAVs('ALL_PSMLT',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
pgaut=getDGAVs('ALL_PGAUT',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
rgacs=getDGAVs('ALL_RGACS',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
rsacr=getDGAVs('ALL_RSACR',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
rsacs=getDGAVs('ALL_RSACS',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);

n=getDGAVs('ALL_Q09',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
q=getDGAVS('ALL_Q04',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);

Rmass_S=n./q;

dq = riacr_s+rsaut+rsbrk - (pssub+psmlt+pgaut).*Rmass_S - rgacs - rsacr - rsacs;

dq2=getDGAVs('ALL_DQ09',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);

case 'snow'
psaut=getDGAVs('ALL_PSAUT',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
psdep=getDGAVs('ALL_PSDEP',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
psacw=getDGAVs('ALL_PSACW',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
psaci=getDGAVs('ALL_PSACI',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
praci_s=getDGAVs('ALL_PRACI_S',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
piacr_s=getDGAVs('ALL_PIACR_S',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
pssub=getDGAVs('ALL_PSSUB',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
pgacs=getDGAVs('ALL_PGACS',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
pracs=getDGAVs('ALL_PRACS',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
pgaut=getDGAVs('ALL_PGAUT',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);
psmlt=getDGAVs('ALL_PSMLT',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);

dgs={'psaut','psdep','psacw','psaci','praci_s','piacr_s',...
    'pssub','pgacs','pracs','pgaut','psmlt'};


dq = psaut+psdep+psacw+psaci+praci_s+piacr_s - ...
    (pssub+pgacs+pracs+pgaut+psmlt);

dq2=getDGAVs('ALL_DQ04',dgstrDan(1).dg,TimeAvDan(1).DGAV,1,1);

end

%ii=47;
%fprintf(1,'%e\n = %e\n = %e\n = %e\n = %e\n = %e\n = %e\n = %e\n = %e\n = %e\n = %e\n',psaut(ii),psdep(ii),psacw(ii),psaci(ii),praci_s(ii),piacr_s(ii), ...
%    pssub(ii),pgacs(ii),pracs(ii),pgaut(ii),psmlt(ii));


