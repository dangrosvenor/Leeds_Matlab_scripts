%clear dgs;

 if  jj==fnmin %only search for column numbers on first pass       
    
    
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
        
        
        for iprc=1:length(dgs);
            nam=strcat('ALL_',dgs{iprc});
            dgcol2(j).d(iprc)=getDGcol(nam,dgstrDan(jc).dg);
            %icediag(j).i(:,jj,iprc)=getDGAVs(nam,dgstrDan(jc).dg); 
        end
        

    
     id=1;        

    dgcol(j).d(id)=getDGcol('ALL_WQ01',dgstrDan(jc).dg);   %getDGcol doesn't divide by the area of the columns
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQ02',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQ03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQ04',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQ05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQ06',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQ07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQ08',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQ09',dgstrDan(jc).dg);
    id=id+1;
    %9
    
    dgcol(j).d(id)=getDGcol('ALL_WQSG01',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQSG02',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQSG03',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQSG04',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQSG05',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQSG06',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQSG07',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQSG08',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_WQSG09',dgstrDan(jc).dg); 
    id=id+1;
    %18
    
      
    dgcol(j).d(id)=getDGcol('ALL_FQ01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_FQ02',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALL_FQ03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_FQ04',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALL_FQ05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_FQ06',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALL_FQ07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_FQ08',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALL_FQ09',dgstrDan(jc).dg);
    id=id+1;
    %27

    
    dgcol(j).d(id)=getDGcol('ALL_DQ01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_DQ02',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_DQ03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_DQ04',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_DQ05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_DQ06',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_DQ07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_DQ08',dgstrDan(jc).dg);
    id=id+1;  
    dgcol(j).d(id)=getDGcol('ALL_DQ09',dgstrDan(jc).dg);
    id=id+1; 
    %36
    
    dgcol(j).d(id)=getDGcol('ALL_Q01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_Q02',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_Q03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_Q04',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_Q05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_Q06',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_Q07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALL_Q08',dgstrDan(jc).dg);
    id=id+1;  
    dgcol(j).d(id)=getDGcol('ALL_Q09',dgstrDan(jc).dg);
    id=id+1;
    %45
    
    %up ...........................
    
    dgcol(j).d(id)=getDGcol('ALu_WQ01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQ02',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQ03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQ04',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQ05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQ06',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQ07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQ08',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQ09',dgstrDan(jc).dg);
    id=id+1;
    %54
    
    dgcol(j).d(id)=getDGcol('ALu_WQSG01',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQSG02',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQSG03',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQSG04',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQSG05',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQSG06',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQSG07',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQSG08',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_WQSG09',dgstrDan(jc).dg); 
    id=id+1;
    %63
      
    dgcol(j).d(id)=getDGcol('ALu_FQ01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_FQ02',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALu_FQ03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_FQ04',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALu_FQ05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_FQ06',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALu_FQ07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_FQ08',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALu_FQ09',dgstrDan(jc).dg);
    id=id+1;
    %72

    
    dgcol(j).d(id)=getDGcol('ALu_DQ01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_DQ02',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_DQ03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_DQ04',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_DQ05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_DQ06',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_DQ07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_DQ08',dgstrDan(jc).dg);
    id=id+1;  
    dgcol(j).d(id)=getDGcol('ALu_DQ09',dgstrDan(jc).dg);
    id=id+1; 
    %81
    
    dgcol(j).d(id)=getDGcol('ALu_Q01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_Q02',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_Q03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_Q04',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_Q05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_Q06',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_Q07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALu_Q08',dgstrDan(jc).dg);
    id=id+1;  
    dgcol(j).d(id)=getDGcol('ALu_Q09',dgstrDan(jc).dg);
    id=id+1;
    %90
    
        %down..........................
        
    
    dgcol(j).d(id)=getDGcol('ALd_WQ01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQ02',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQ03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQ04',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQ05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQ06',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQ07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQ08',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQ09',dgstrDan(jc).dg);
    id=id+1;
    %99
    
    dgcol(j).d(id)=getDGcol('ALd_WQSG01',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQSG02',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQSG03',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQSG04',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQSG05',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQSG06',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQSG07',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQSG08',dgstrDan(jc).dg); 
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_WQSG09',dgstrDan(jc).dg); 
    id=id+1;
    %108
      
    dgcol(j).d(id)=getDGcol('ALd_FQ01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_FQ02',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALd_FQ03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_FQ04',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALd_FQ05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_FQ06',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALd_FQ07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_FQ08',dgstrDan(jc).dg);
    id=id+1;
     dgcol(j).d(id)=getDGcol('ALd_FQ09',dgstrDan(jc).dg);
    id=id+1;
    %117

    
    dgcol(j).d(id)=getDGcol('ALd_DQ01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_DQ02',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_DQ03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_DQ04',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_DQ05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_DQ06',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_DQ07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_DQ08',dgstrDan(jc).dg);
    id=id+1;  
    dgcol(j).d(id)=getDGcol('ALd_DQ09',dgstrDan(jc).dg);
    id=id+1; 
    %126
    
    dgcol(j).d(id)=getDGcol('ALd_Q01',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_Q02',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_Q03',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_Q04',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_Q05',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_Q06',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_Q07',dgstrDan(jc).dg);
    id=id+1;
    dgcol(j).d(id)=getDGcol('ALd_Q08',dgstrDan(jc).dg);
    id=id+1;  
    dgcol(j).d(id)=getDGcol('ALd_Q09',dgstrDan(jc).dg);
    id=id+1;
    %135
    
end
% 
    for idg=1:length(dgcol(j).d)
        icediagALL(j).i(:,jj,idg)=TimeAvDan(j).DGAV(:,dgcol(j).d(idg));
    end
    
     for idg=1:length(dgcol2(j).d)
         icediag(j).i(:,jj,idg)=TimeAvDan(j).DGAV(:,dgcol2(j).d(idg));
     end
    
    
%     w_maxprof(j).w(:,jj)=max(TwoDDan(j).W,[],2);
%         w_minprof(j).w(:,jj)=min(TwoDDan(j).W,[],2);

%         prcs=[0:5:100];
%         vap_prctiles(j).v(:,jj,1:length(prcs))=(prctile(TwoDDan(j).Q(:,:,1)',prcs))';
%         time(j,jj)=SerDan(j).SER(end,1);