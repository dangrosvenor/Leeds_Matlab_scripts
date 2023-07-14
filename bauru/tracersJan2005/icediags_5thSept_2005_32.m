%clear dgs;


i3d=1;

microcase='oldnew'; %csae for original process rates
microcase='numrates'; %case for the new process rates with the number rates scaled to match the acutal microphysical change


if  jjfile==fnmin %only search for column numbers on first pass       
%     
    clear dgs dgcol2
    
        switch microcase
        case 'oldold'
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
			dgs{16}='PSDEP';
			dgs{17}='PIACR_S';
			dgs{18}='PSACW';
			dgs{19}='PSSUB';
			dgs{20}='PGACS';
			dgs{21}='PRACS';
			dgs{22}='PGAUT';
			dgs{23}='PSMLT';
			dgs{24}='RSBRK';
			dgs{25}='RIACR_S';
			dgs{26}='RGACS';
            dgs{27}='RSACR';
            dgs{28}='RSACS';
            
            dgs{29}='PGSUB';
            
        case 'oldnew'
            dgs={'PGDEP','PGMLT', ...  %1-2
          'PRAUT',   'PGSHD',   'PRACW',   'PSMLT', ... %3-6
          'PIMLT',   'PSAUT',   'PSDEP',   'PIACR_G', ... %7-10
          'PSACI',   'PGACI',   'PGACW',   'PGACS', ... %11-14
          'PGACR',   'PSACR',   'PRACS',   'PSACW', ... %15-18
          'PGAUT',   'PGFR',    'PRACI_G', 'PGWET', ... %19-22
          'PGDRY',   'PGSUB',   'PSSUB',   'PREVP', ... %23-26
          'PISUB',   'DQI  ',   'PIHAL',   'PIPRM', ... %27-30
          'PIDEP',   'PIACW',   'PICNT',   'PIFRW', ... %31-34
          'PIACR_S', 'PRACI_S', 'PRACI',   'PIACR', ... %35-38
          'RIACR_S', 'RSACR',   'RIACR_G', 'RGFR', ...  %39-42
          'RGACS',   'RGACR',   'RSAUT',   'RIACI', ... %43-46
          'RSACS',   'RSBRK'  ...                       %47-48
          };
            
            imrsour=[29:34]; %sources of ice mixing ratio 
            imrsink=[8 11 12 27 7 21 36];
            imrnuc=[29 30 33 34];
            imrnonmpc=[11 12 21 36 29 32]; %accretion type processes not included in the MPC
            imrmpc=[8 27 7 30 31 33 34]; %accretion type processes not included in the MPC
            
            smrsour=[8 9 18 11 35 36]; %sources of snow mixing ratio 
            smrsink=[25 14 17 19 6];
            smrmpc=[8 9 25 19 6]; %processes included in the MPC
            
            gmrsour=[1 12:17 10 19:21]; %sources of graupel mixing ratio 
            gmrsink=[2 24 4];
            
            hmssour=[1 13 15 20 16 10 imrsour 9 18 35]; %sources of all ice species - i.e. not including ones that convert from one ice to another
            hmssink=[2 4 6 7 24 25 27]; %as above but sinks
            
            liqsink=[13 3 18 5 32 29 33 34];
            
            
        case 'numrates'
            dgs={'PGDEP','PGMLT', ...  %1-2
          'PRAUT',   'PGSHD',   'PRACW',   'PSMLT', ... %3-6
          'PIMLT',   'PSAUT',   'PSDEP',   'PIACR_G', ... %7-10
          'PSACI',   'PGACI',   'PGACW',   'PGACS', ... %11-14
          'PGACR',   'PSACR',   'PRACS',   'PSACW', ... %15-18
          'PGAUT',   'PGFR',    'PRACI_G', 'PGWET', ... %19-22
          'PGDRY',   'PGSUB',   'PSSUB',   'PREVP', ... %23-26
          'PISUB',   'DQI  ',   'PIHAL',   'PIPRM', ... %27-30
          'PIDEP',   'PIACW',   'PICNT',   'PIFRW', ... %31-34
          'PIACR_S', 'PRACI_S', 'PRACI',   'PIACR', ... %35-38
          'RIACR_S', 'RSACR',   'RIACR_G', 'RGFR', ...  %39-42
          'RGACS',   'RGACR',   'RSAUT',   'RIACI', ... %43-46
          'RSACS',   'RSBRK'  ...                       %47-48
          , 'RGmlt'		  ... %49
          , 'RSmlt'		   ...%50
          , 'RImlt'		  ...
          , 'RSaci'		  ...
          , 'RGaci'		  ... 
          , 'RGacw'		  ...  %??? - not used?
          , 'RRacs'		  ... %55
          , 'RGaut'		  ...
          , 'RGsub'		  ... 
          , 'RSsub'		  ...
          , 'RIcnt'		  ...
          , 'RIhal'		  ... %60
          , 'RIprm'		  ... 
          , 'RIfrw'		  ... 
          , 'RRaci_S'	  ... 
          , 'RRaci_G'	  ...
          , 'RSaut_I'	  ... %65
          , 'RSaut_S'	  ...
          , 'RGaut_S'	  ... 
          , 'RGaut_G'	  ...
          , 'RSacr_S'	  ...
          , 'RSacr_G' }   ;   %70
        end
        
        ivapsource=[24:27];
        ivapsink=  [9 1 31 30]; %for oldnew and numrates options 
        
        icencsour=[59:62]; incmpc=[59 61 62 45]; %PIPRM - heterogeneous nucleation not in MPC, but homog would replace it
        icencsink=[52 53 63 64 45 46];
        
        snowncsour=[39 45 48]; sncmpc=[45 48 58 50 56]; %not sure whether RSBRK(48) is in MPC
        snowncsink=[58 50 56 43 40 47];
        
        gncsour=[56 40 41 42];
        gncsink=[57 49];
        
        error=0;
        for iprc=1:length(dgs);
            nam=strcat('ALL_',dgs{iprc});
            [dgcol2(j).d(iprc),error]=getDGcol(nam,dgstrDan(jc).dg,error);
            %icediag(j).i(:,jj,iprc)=getDGAVs(nam,dgstrDan(jc).dg,error); 
        end
        
        if error>=1
            fprintf(1,'\n*********** WARNING %d calls to getDGcols for microphys failed *************\n',error);
        end
        
% 
    error=0;

     id=1;        

    [dgcol(j).d(id),error]=getDGcol('ALL_WQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ09',dgstrDan(jc).dg,error);
    id=id+1;
    %10 %next one will be number 10
    
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG01',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG02',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG03',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG04',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG05',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG06',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG07',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG08',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG09',dgstrDan(jc).dg,error); 
    id=id+1;
    %19
    
      
    [dgcol(j).d(id),error]=getDGcol('ALL_FQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_FQ02',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALL_FQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_FQ04',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALL_FQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_FQ06',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALL_FQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_FQ08',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALL_FQ09',dgstrDan(jc).dg,error);
    id=id+1;
    %28

    
    [dgcol(j).d(id),error]=getDGcol('ALL_DQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_DQ02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_DQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_DQ04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_DQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_DQ06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_DQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_DQ08',dgstrDan(jc).dg,error);
    id=id+1;  
    [dgcol(j).d(id),error]=getDGcol('ALL_DQ09',dgstrDan(jc).dg,error);
    id=id+1; 
    %37
    
    [dgcol(j).d(id),error]=getDGcol('ALL_Q01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q08',dgstrDan(jc).dg,error);
    id=id+1;  
    [dgcol(j).d(id),error]=getDGcol('ALL_Q09',dgstrDan(jc).dg,error);
    id=id+1;
    %46
    
    %up ...........................
    
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ09',dgstrDan(jc).dg,error);
    id=id+1;
    %55
    
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG01',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG02',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG03',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG04',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG05',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG06',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG07',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG08',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG09',dgstrDan(jc).dg,error); 
    id=id+1;
    %64
      
    [dgcol(j).d(id),error]=getDGcol('ALu_FQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_FQ02',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALu_FQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_FQ04',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALu_FQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_FQ06',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALu_FQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_FQ08',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALu_FQ09',dgstrDan(jc).dg,error);
    id=id+1;
    %73

    
    [dgcol(j).d(id),error]=getDGcol('ALu_DQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_DQ02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_DQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_DQ04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_DQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_DQ06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_DQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_DQ08',dgstrDan(jc).dg,error);
    id=id+1;  
    [dgcol(j).d(id),error]=getDGcol('ALu_DQ09',dgstrDan(jc).dg,error);
    id=id+1; 
    %82
    
    [dgcol(j).d(id),error]=getDGcol('ALu_Q01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q08',dgstrDan(jc).dg,error);
    id=id+1;  
    [dgcol(j).d(id),error]=getDGcol('ALu_Q09',dgstrDan(jc).dg,error);
    id=id+1;
    %91
    
        %down..........................
        
    
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ09',dgstrDan(jc).dg,error);
    id=id+1;
    %100
    
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG01',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG02',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG03',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG04',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG05',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG06',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG07',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG08',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG09',dgstrDan(jc).dg,error); 
    id=id+1;
    %109
      
    [dgcol(j).d(id),error]=getDGcol('ALd_FQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_FQ02',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALd_FQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_FQ04',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALd_FQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_FQ06',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALd_FQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_FQ08',dgstrDan(jc).dg,error);
    id=id+1;
     [dgcol(j).d(id),error]=getDGcol('ALd_FQ09',dgstrDan(jc).dg,error);
    id=id+1;
    %118

    
    [dgcol(j).d(id),error]=getDGcol('ALd_DQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_DQ02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_DQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_DQ04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_DQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_DQ06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_DQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_DQ08',dgstrDan(jc).dg,error);
    id=id+1;  
    [dgcol(j).d(id),error]=getDGcol('ALd_DQ09',dgstrDan(jc).dg,error);
    id=id+1; 
    %127
    
    [dgcol(j).d(id),error]=getDGcol('ALd_Q01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q08',dgstrDan(jc).dg,error);
    id=id+1;  
    [dgcol(j).d(id),error]=getDGcol('ALd_Q09',dgstrDan(jc).dg,error);
    id=id+1;
    %136
    
%new ones added after Nov 2005

    [dgcol(j).d(id),error]=getDGcol('ALL_W',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_W',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_W',dgstrDan(jc).dg,error);
    id=id+1;
    
    %139
    
%below added in Jan 2006 to get tracer values

    [dgcol(j).d(id),error]=getDGcol('ALL_WQ10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQ15',dgstrDan(jc).dg,error);
    id=id+1;
    
    %145
    
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG10',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG11',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG12',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG13',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG14',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_WQSG15',dgstrDan(jc).dg,error); 
    id=id+1;
  
    %151
    
    [dgcol(j).d(id),error]=getDGcol('ALL_Q10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q11',dgstrDan(jc).dg,error);
    id=id+1;   
    [dgcol(j).d(id),error]=getDGcol('ALL_Q12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALL_Q15',dgstrDan(jc).dg,error);
    id=id+1;
  
    %157
    
    %up ...........................
    
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQ15',dgstrDan(jc).dg,error);
    id=id+1;
  
    %163
    
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG10',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG11',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG12',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG13',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG14',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_WQSG15',dgstrDan(jc).dg,error); 
    id=id+1;
 
    %169

    
    [dgcol(j).d(id),error]=getDGcol('ALu_Q10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALu_Q15',dgstrDan(jc).dg,error);
    id=id+1;
   
     %175;
    
        %down..........................
        
    
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQ15',dgstrDan(jc).dg,error);
    id=id+1;
  
	%181
    
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG10',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG11',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG12',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG13',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG14',dgstrDan(jc).dg,error); 
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_WQSG15',dgstrDan(jc).dg,error); 
    id=id+1;
 
    %185
    
    [dgcol(j).d(id),error]=getDGcol('ALd_Q10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_Q15',dgstrDan(jc).dg,error);
    id=id+1;
    
    %193
    
% ACC - all cloudy points ****************************************

    
   [dgcol(j).d(id),error]=getDGcol('ACC_WQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ09',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQ15',dgstrDan(jc).dg,error);
    id=id+1;
    
    %208
    
    
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG09',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_WQSG15',dgstrDan(jc).dg,error);
    id=id+1;
    
   %223
 
    
    [dgcol(j).d(id),error]=getDGcol('ACC_Q01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q09',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_Q15',dgstrDan(jc).dg,error);
    id=id+1;  
    %238
    
    [dgcol(j).d(id),error]=getDGcol('ALL_LW',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ALL_SW',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ACC_LW',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ACC_SW',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ALu_LW',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ALu_SW',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ALd_LW',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ALd_SW',dgstrDan(jc).dg,error);
	id=id+1;
    %246
    
	[dgcol(j).d(id),error]=getDGcol('ALL_TH',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ACC_TH',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ALu_TH',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ALd_TH',dgstrDan(jc).dg,error);
	id=id+1;
    %250 next
    
    [dgcol(j).d(id),error]=getDGcol('CLu_Q01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q09',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_Q15',dgstrDan(jc).dg,error);
    id=id+1;
    %265 next (add on no. just found dgcol for
    
    [dgcol(j).d(id),error]=getDGcol('W>1_Q01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q09',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_Q15',dgstrDan(jc).dg,error);
    id=id+1;
	%280 next (add on no. just found dgcol for
    
    [dgcol(j).d(id),error]=getDGcol('ALu_A',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ALd_A',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLu_A',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_A',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACC_A',dgstrDan(jc).dg,error);
    id=id+1; 
    [dgcol(j).d(id),error]=getDGcol('ACu_A',dgstrDan(jc).dg,error);
    id=id+1; 
	%286 next (add on no. just found dgcol for
    
    [dgcol(j).d(id),error]=getDGcol('ACu_Q01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q09',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_Q15',dgstrDan(jc).dg,error);
    id=id+1;
    %301 next (add on no. just found dgcol for
    
    %new ones added after Sep 2006

    [dgcol(j).d(id),error]=getDGcol('CLu_W',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_W',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('W>1_W',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('CLd_W',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_W',dgstrDan(jc).dg,error);
    id=id+1;
    
    
    
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ09',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACu_WQ15',dgstrDan(jc).dg,error);
    id=id+1;
    
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ09',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_WQ15',dgstrDan(jc).dg,error);
    id=id+1;
    
    [dgcol(j).d(id),error]=getDGcol('ACd_Q01',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q02',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q03',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q04',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q05',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q06',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q07',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q08',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q09',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q10',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q11',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q12',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q13',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q14',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('ACd_Q15',dgstrDan(jc).dg,error);
    id=id+1; 
    
    [dgcol(j).d(id),error]=getDGcol('ACd_A',dgstrDan(jc).dg,error);
    id=id+1; 
    [dgcol(j).d(id),error]=getDGcol('CLd_A',dgstrDan(jc).dg,error);
    id=id+1;
    
     %next one = 353  (add on no. just found dgcol for)
        
    [dgcol(j).d(id),error]=getDGcol('WTHSG',dgstrDan(jc).dg,error);
    id=id+1; 
    [dgcol(j).d(id),error]=getDGcol('WTHAD',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('VW',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('VWSG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WKE',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WKESG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('VV',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('VVSG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WW',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WWSG',dgstrDan(jc).dg,error);
    id=id+1;
    
    %next one = 363  (add on no. just found dgcol for)
    
    [dgcol(j).d(id),error]=getDGcol('WQ01AD',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ02AD',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ03AD',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ04AD',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ05AD',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ06AD',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ07AD',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ08AD',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ09AD',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ10AD',dgstrDan(jc).dg,error);
    id=id+1;
    
    %next one = 373  (add on no. just found dgcol for)
    
    [dgcol(j).d(id),error]=getDGcol('WQ01SG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ02SG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ03SG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ04SG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ05SG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ06SG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ07SG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ08SG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ09SG',dgstrDan(jc).dg,error);
    id=id+1;
    [dgcol(j).d(id),error]=getDGcol('WQ10SG',dgstrDan(jc).dg,error);
    
    %next one = 383  (add on no. just found dgcol for)    
    
	[dgcol(j).d(id),error]=getDGcol('ALL_TEMP',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ACC_TEMP',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ALu_TEMP',dgstrDan(jc).dg,error);
	id=id+1;
	[dgcol(j).d(id),error]=getDGcol('ALd_TEMP',dgstrDan(jc).dg,error);
	id=id+1;
    %250 next

    




    
    
    if error>=1
        fprintf(1,'\n*********** WARNING %d calls to getDGcols for icediagsALL failed *************\n',error);
    end
    
end
% 
    for idg=1:length(dgcol(j).d)
        icediagsALL(j).i(:,jj,idg)=TimeAvDan(jc).DGAV(:,dgcol(j).d(idg)); %saved as icediagsALL_a-b from Nov 2005
    end
    
     for idg=1:length(dgcol2(j).d)
%          switch microcase
%          case 'old'
             icediag(j).i(:,jj,idg)=TimeAvDan(jc).DGAV(:,dgcol2(j).d(idg));
%          case 'oldnew'
%              icediag(j).i(:,jj,idg)=TimeAvDan(jc).DGAV(:,dgcol2(j).d(idg)); %saved as icediag_a-b from Nov 2005
%          case 'numrates'
%              icediags_nums(j).i(:,jj,idg)=TimeAvDan(jc).DGAV(:,dgcol2(j).d(idg));
%          end             
     end
    
    if i3d~=1
         MaxW(j).w(:,jj)=max(TwoD.W,[],2);
         MinW(j).w(:,jj)=min(TwoD.W,[],2);
         
         prc=[0:5:100];
%         w_prctiles(j).w(1:length(prc),1:length(GridDan(jc).Z),jj)=prctile(TwoD.W',prc);
         w_prctiles(j).w(1:length(prc),1:length(GridDan(jc).Z),jj)=prctile(TwoD.W',prc);
         
      for iq=1:size(TwoD.Q,3)   
         q_maxprof(j).q(:,jj,iq)=max(TwoD.Q(:,:,iq),[],2);
         q_prctiles(j).q(1:length(prc),1:length(GridDan(jc).Z),jj,iq)=prctile(TwoD.Q(:,:,iq)',prc);
      end
      
                 
         
         
         
         
         
         if jjfile==fnmin %set up extreme max and min on first timestep
             minov=9999;
             maxov=-9999;
             clear wbb
         end
         
         maxov=ceil(max([maxov maxALL(w_prctiles(j).w(:,:,jj))])); %update max and min
         minov=floor(min([minov minALL(w_prctiles(j).w(:,:,jj))]));
         
         step=0.5;
         wb=[minov:step:maxov];
		%dw=wb(2:end)-wb(1:end-1);
        ww=TwoD.W;		
        
		ntot=size(ww,2)*size(ww,3); %total number of points in population, each is 5% of that point in terms of number
		
		for i=2:length(wb)
            a=find(ww>=wb(i-1) & ww<wb(i));
            wdist(j).w(i-1,jj)=length(a)/step/ntot; %total number of points in 2d slice in this bin per dw and normalised for all points in slice
		end
        wbb(j).w(1:length(wb),jj)=wb; %store w bins as they will vary

%         prcs=[0:5:100];
%         vap_prctiles(j).v(:,jj,1:length(prcs))=(prctile(TwoDDan(j).Q(:,:,1)',prcs))';
%         time(j,jj)=SerDan(j).SER(end,1);

			tref=repmat(GridDan(idir).THREF,[1 length(GridDan(idir).Y1)]);
            T=squeeze(sum(TwoD.TH1,3))+tref;
            P=TwoD.PP;
            T=T./(1e5./P).^0.286;
            qsi=satvapPress(T,'buck2','ice',P,1)/f; %satvappress gives in ppmv if final flag=1
            si=100*(TwoD.Q(:,:,1)-qsi)./qsi;
            
            simaxTimH(1).s(:,jj)=max(si,[],2);
            siminTimH(1).s(:,jj)=min(si,[],2);
            simean(1).s(:,jj)=mean(si,2);
            
            
      end   %i3d
