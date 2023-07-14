clear dgs;

 if  jj==fnmin %only search for column numbers on first pass       
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
    id=1;        
 
%     dgs2{id}='ALu_WQ01';
%     id=id+1;
%     dgs2{id}='ALd_WQ01';
%     id=id+1;
%     dgs2{id}='ALu_WQSG01'; 
%     id=id+1;
%     dgs2{id}='ALd_WQSG01';
%     id=id+1;
%     dgs2{id}='ALu_FQ01';
%     id=id+1;
%     dgs2{id}='ALd_FQ01';
%     id=id+1;
%     
%     
%     dgs2{id}='ALu_WQ04';
%     id=id+1;
%     dgs2{id}='ALd_WQ04';
%     id=id+1;
%     dgs2{id}='ALu_WQSG04'; 
%     id=id+1;
%     dgs2{id}='ALd_WQSG04';
%     id=id+1;
%     dgs2{id}='ALu_FQ04';
%     id=id+1;
%     dgs2{id}='ALd_FQ04';
%     id=id+1;
%     
%     dgs2{id}='ALu_WQ05';
%     id=id+1;
%     dgs2{id}='ALd_WQ05';
%     id=id+1;
%     dgs2{id}='ALu_WQSG05'; 
%     id=id+1;
%     dgs2{id}='ALd_WQSG05';
%     id=id+1;
%     dgs2{id}='ALu_FQ05';
%     id=id+1;
%     dgs2{id}='ALd_FQ05';
%     id=id+1;
%     
%     dgs2{id}='ALu_WQ06';
%     id=id+1;
%     dgs2{id}='ALd_WQ06';
%     id=id+1;
%     dgs2{id}='ALu_WQSG06'; 
%     id=id+1;
%     dgs2{id}='ALd_WQSG06';
%     id=id+1;
%     dgs2{id}='ALu_FQ06';
%     id=id+1;
%     dgs2{id}='ALd_FQ06';
%     id=id+1;
%     
%     dgs2{id}='ALu_A';
%     id=id+1;
%     dgs2{id}='ALd_A';
%     id=id+1;
%     
%     dgs2{id}='ALu_DQ01';
%     id=id+1;
%     dgs2{id}='ALd_DQ01';
%     id=id+1;
%     dgs2{id}='ALu_DQ04';
%     id=id+1;
%     dgs2{id}='ALd_DQ04';
%     id=id+1;
%     dgs2{id}='ALu_DQ05';
%     id=id+1;
%     dgs2{id}='ALd_DQ05';
%     id=id+1;
%     dgs2{id}='ALu_DQ06';
%     id=id+1;
%     dgs2{id}='ALd_DQ06';
%     id=id+1;  
%     
%     dgs2{id}='ALu_Q01';
%     id=id+1;
%     dgs2{id}='ALd_Q01';
%     id=id+1;
%     dgs2{id}='ALu_Q04';
%     id=id+1;
%     dgs2{id}='ALd_Q04';
%     id=id+1;
%     dgs2{id}='ALu_Q05';
%     id=id+1;
%     dgs2{id}='ALd_Q05';
%     id=id+1;
%     dgs2{id}='ALu_Q06';
%     id=id+1;
%     dgs2{id}='ALd_Q06';
%     id=id+1;
    %%%%
    
    dgs2{id}='ALu_WQ07';
    id=id+1;
    dgs2{id}='ALd_WQ07';
    id=id+1;
    dgs2{id}='ALu_WQSG07'; 
    id=id+1;
    dgs2{id}='ALd_WQSG07';
    id=id+1;
    dgs2{id}='ALu_FQ07';
    id=id+1;
    dgs2{id}='ALd_FQ07';
    id=id+1;
    
    
    dgs2{id}='ALu_WQ08';
    id=id+1;
    dgs2{id}='ALd_WQ08';
    id=id+1;
    dgs2{id}='ALu_WQSG08'; 
    id=id+1;
    dgs2{id}='ALd_WQSG08';
    id=id+1;
    dgs2{id}='ALu_FQ08';
    id=id+1;
    dgs2{id}='ALd_FQ08';
    id=id+1;
    
    dgs2{id}='ALu_WQ09';
    id=id+1;
    dgs2{id}='ALd_WQ09';
    id=id+1;
    dgs2{id}='ALu_WQSG09'; 
    id=id+1;
    dgs2{id}='ALd_WQSG09';
    id=id+1;
    dgs2{id}='ALu_FQ09';
    id=id+1;
    dgs2{id}='ALd_FQ09';
    id=id+1;
    
    dgs2{id}='ALu_WQ02';
    id=id+1;
    dgs2{id}='ALd_WQ02';
    id=id+1;
    dgs2{id}='ALu_WQSG02'; 
    id=id+1;
    dgs2{id}='ALd_WQSG02';
    id=id+1;
    dgs2{id}='ALu_FQ02';
    id=id+1;
    dgs2{id}='ALd_FQ02';
    id=id+1;
    
    
    dgs2{id}='ALu_WQ03';
    id=id+1;
    dgs2{id}='ALd_WQ03';
    id=id+1;
    dgs2{id}='ALu_WQSG03'; 
    id=id+1;
    dgs2{id}='ALd_WQSG03';
    id=id+1;
    dgs2{id}='ALu_FQ03';
    id=id+1;
    dgs2{id}='ALd_FQ03';
    id=id+1;
    
 
    
    dgs2{id}='ALu_DQ07';
    id=id+1;
    dgs2{id}='ALd_DQ07';
    id=id+1;
    dgs2{id}='ALu_DQ08';
    id=id+1;
    dgs2{id}='ALd_DQ08';
    id=id+1;
    dgs2{id}='ALu_DQ09';
    id=id+1;
    dgs2{id}='ALd_DQ09';
    id=id+1;
    dgs2{id}='ALu_DQ02';
    id=id+1;
    dgs2{id}='ALd_DQ02';
    id=id+1;  
    dgs2{id}='ALu_DQ03';
    id=id+1;
    dgs2{id}='ALd_DQ03';
    id=id+1;  
    
    dgs2{id}='ALu_Q07';
    id=id+1;
    dgs2{id}='ALd_Q07';
    id=id+1;
    dgs2{id}='ALu_Q08';
    id=id+1;
    dgs2{id}='ALd_Q08';
    id=id+1;
    dgs2{id}='ALu_Q09';
    id=id+1;
    dgs2{id}='ALd_Q09';
    id=id+1;
    dgs2{id}='ALu_Q02';
    id=id+1;
    dgs2{id}='ALd_Q02';
    id=id+1;
    dgs2{id}='ALu_Q03';
    id=id+1;
    dgs2{id}='ALd_Q03';
    id=id+1;
    
    for iprc=1:length(dgs2);
        nam=dgs2{iprc};
        dgcol(j).d(iprc)=getDGcol(nam,dgstrDan(jc).dg);
        %icediag(j).i(:,jj,iprc)=getDGAVs(nam,dgstrDan(jc).dg); 
    end

end
% 
     for idg=1:length(dgcol(j).d)
         icediag5(j).i(:,jj,idg)=TimeAvDan(j).DGAV(:,dgcol(j).d(idg));
     end
%     
%     for idg=1:length(dgcol2(j).d)
%         icediag(j).i(:,jj,idg)=TimeAvDan(j).DGAV(:,dgcol2(j).d(idg));
%     end
    
    
%    w_maxprof(j).w(:,jj)=max(TwoDDan(j).W,[],2);
%        w_minprof(j).w(:,jj)=min(TwoDDan(j).W,[],2);

%         prcs=[0:5:100];
%         vap_prctiles(j).v(:,jj,1:length(prcs))=(prctile(TwoDDan(j).Q(:,:,1)',prcs))';
%         time(j,jj)=SerDan(j).SER(end,1);