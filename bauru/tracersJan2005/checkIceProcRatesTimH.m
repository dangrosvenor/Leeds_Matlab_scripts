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

zmin=16e3;
zmax=22e3;

izmin=find(z<=zmin);
izmin=izmin(end);
izmax=find(z>=zmax);
izmax=izmax(1);

tlim=[16.8 16.9];
it1=find(time<=tlim(1));it1=it1(end);
it2=find(time>=tlim(2));it2=it2(1);

for i=1:12
    ma=max(max(icediag(1).i(izmin:izmax,it1:it2,i)));
     fprintf(1,'(%d) Prate %s : max = %e\n',i,dgs{i},ma);
 end