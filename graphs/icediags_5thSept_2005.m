clear dgs;

% if  jj==fnmin %only search for column numbers on first pass       
%     
%     
%         dgs{1}='PIMLT';
% 		dgs{2}='PSAUT';
% 		dgs{3}='PSACI';
% 		dgs{4}='PRACI_S';
% 		dgs{5}='PGACI';
% 		dgs{6}='PRACI_G';
% 		dgs{7}='PIHAL';
% 		dgs{8}='PIPRM';
% 		dgs{9}='PICNT';
% 		dgs{10}='PIDEP';
% 		dgs{11}='PIACW';
% 		dgs{12}='PIFRW';
%         dgs{13}='RSAUT';
%         dgs{14}='RIACI';
%         
%         dgs{15}='PISUB';
% 		dgs{16}='PSDEP';
% 		dgs{17}='PIACR_S';
% 		dgs{18}='PSACW';
% 		dgs{19}='PSSUB';
% 		dgs{20}='PGACS';
% 		dgs{21}='PRACS';
% 		dgs{22}='PGAUT';
% 		dgs{23}='PSMLT';
% 		dgs{24}='RSBRK';
% 		dgs{25}='RIACR_S';
% 		dgs{26}='RGACS';
%         dgs{27}='RSACR';
%         dgs{28}='RSACS';
%         
%         dgs{29}='PGSUB';
%         
%         
%         for iprc=1:length(dgs);
%             nam=strcat('ALL_',dgs{iprc});
%             dgcol2(j).d(iprc)=getDGcol(nam,dgstrDan(jc).dg);
%             %icediag(j).i(:,jj,iprc)=getDGAVs(nam,dgstrDan(jc).dg); 
%         end
%         
% 
%     
%     id=1;        
% 
%     dgcol(j).d(id)=getDGcol('ALu_WQ01',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_WQ01',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_WQSG01',dgstrDan(jc).dg); 
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_WQSG01',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_FQ01',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_FQ01',dgstrDan(jc).dg);
%     id=id+1;
%     
%     
%     dgcol(j).d(id)=getDGcol('ALu_WQ04',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_WQ04',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_WQSG04',dgstrDan(jc).dg); 
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_WQSG04',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_FQ04',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_FQ04',dgstrDan(jc).dg);
%     id=id+1;
%     
%     dgcol(j).d(id)=getDGcol('ALu_WQ05',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_WQ05',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_WQSG05',dgstrDan(jc).dg); 
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_WQSG05',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_FQ05',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_FQ05',dgstrDan(jc).dg);
%     id=id+1;
%     
%     dgcol(j).d(id)=getDGcol('ALu_WQ06',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_WQ06',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_WQSG06',dgstrDan(jc).dg); 
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_WQSG06',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_FQ06',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_FQ06',dgstrDan(jc).dg);
%     id=id+1;
%     
%     dgcol(j).d(id)=getDGcol('ALu_A',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_A',dgstrDan(jc).dg);
%     id=id+1;
%     
%     dgcol(j).d(id)=getDGcol('ALu_DQ01',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_DQ01',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_DQ04',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_DQ04',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_DQ05',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_DQ05',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_DQ06',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_DQ06',dgstrDan(jc).dg);
%     id=id+1;  
%     
%     dgcol(j).d(id)=getDGcol('ALu_Q01',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_Q01',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_Q04',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_Q04',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_Q05',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_Q05',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALu_Q06',dgstrDan(jc).dg);
%     id=id+1;
%     dgcol(j).d(id)=getDGcol('ALd_Q06',dgstrDan(jc).dg);
%     id=id+1;  
% end
% 
%     for idg=1:length(dgcol(j).d)
%         icediag4(j).i(:,jj,idg)=TimeAvDan(j).DGAV(:,dgcol(j).d(idg));
%     end
%     
%     for idg=1:length(dgcol2(j).d)
%         icediag(j).i(:,jj,idg)=TimeAvDan(j).DGAV(:,dgcol2(j).d(idg));
%     end
    
    
    w_maxprof(j).w(:,jj)=max(TwoDDan(j).W,[],2);
        w_minprof(j).w(:,jj)=min(TwoDDan(j).W,[],2);

%         prcs=[0:5:100];
%         vap_prctiles(j).v(:,jj,1:length(prcs))=(prctile(TwoDDan(j).Q(:,:,1)',prcs))';
%         time(j,jj)=SerDan(j).SER(end,1);