dgs{1}='PISUB';
dgs{2}='PSDEP';
dgs{3}='PIACR_S';
dgs{4}='PSACW';
dgs{5}='PSSUB';
dgs{6}='PGACS';
dgs{7}='PRACS';
dgs{8}='PGAUT';
dgs{9}='PSMLT';
dgs{10}='RSBRK';
dgs{11}='RIACR_S';
dgs{12}='RGACS';
dgs{13}='RSACR';
dgs{14}='RSACS';
        
exdir='g:/runs/sdlavap/results/processrates/snowprocesses/'

for imp=1:length(dgs)
    titinput=dgs{imp};
    inputted=icediagSnow(1).i(:,:,imp);
    
   plotTimeHeightVap2;
   
   exname=[exdir dgs{imp}];
   set(gcf,'paperpositionmode','auto');
   print(gcf,'-dmeta',exname);
   close(gcf);
   
end