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
dgs{13}='RSAUT';
dgs{14}='RIACI';
dgs{15}='PISUB';
        
exdir='g:/runs/sdlavap/results/processrates/test/';

inorm=[13 14]; %no factor needed
for imp=inorm
    titinput=dgs{imp};
    inputted=icediag(1).i(:,:,imp);
    
   plotTimeHeightVap2;
   
   exname=[exdir dgs{imp}];
   set(gcf,'paperpositionmode','auto');
   print(gcf,'-dmeta',exname);
   close(gcf);
   
end



imnq=[3 4 5 6]; %to be multiplied by n/q
for imp=imnq
    titinput=dgs{imp};
    inputted=icediag(1).i(:,:,imp).*ndivqav;
    
   plotTimeHeightVap2;
   
   exname=[exdir dgs{imp}];
   set(gcf,'paperpositionmode','auto');
   print(gcf,'-dmeta',exname);
   close(gcf);
   
end

MI0=1e-15;
imi0=[7 8 9 12]; %indices for processes to be multiplied by MI0
for imp=imi0
    titinput=dgs{imp};
    inputted=icediag(1).i(:,:,imp)/MI0;
    
   plotTimeHeightVap2;
   
   exname=[exdir dgs{imp}];
   set(gcf,'paperpositionmode','auto');
   print(gcf,'-dmeta',exname);
   close(gcf);
   
end

