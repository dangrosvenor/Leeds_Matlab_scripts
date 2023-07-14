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

zmax=15e3;

zz=find(Grid.Z>zmax);
idir=1;
for i=1:length(dgs)
    dgfind=findhead(dgs{i},dgstrDan(idir).dg);
    dgsf(i)=dgfind(1);
    maxdg=max(TimeAvDan(idir).DGAV(zz(1):end,dgfind(1)));
    fprintf(1,'(%d) Prate %s : max = %e\n',i,dgstrDan(idir).dg{dgfind(1)},maxdg);
end