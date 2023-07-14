function icol=iemm_mphys_cols(dgname,flowing,iflag)
%flowing = 0 : debris region (uses w_deb)
%flowing = 1 : main updraught (uses w_areamean_channel)
%flowing = 2 : outflow region(?) (uses w_out) (scr)
%flowing = 3 : downdraught?

dgs={'ihomog','ievap_num','icond','ievap_size','iauto','ievap_zbase','ievap_ctop','iauto_res','iauto_purge'};
dgs={'ihomog','ievap_num','icond','ievap_size','iauto','ievap_zbase','ievap_ctop','iauto_res','iauto_purge','iprimary','isecondary','ihomog_rain',...
        'ihomog_aero'};

%last one was number 13

if iflag==0

for i=1:length(dgs)
    if strcmp(dgs{i},dgname)==1
        idgs=i;
        break
    end
end

%icol=1 + length(dgs)*2*flowing + idgs*2;
icol =  2 + length(dgs)*flowing + idgs;

else
    
    iflowing=ceil((dgname-2)/(2*length(dgs))) - 1;   
    
    if mod(dgname,2)==0
        icol{3}='nums';
        dgname=dgname-1;
    else
        icol{3}='mass';
    end
    
         
    icol{1}=iflowing;
    idgs=(dgname - 1 - length(dgs)*2*iflowing)/2;
    icol{2}=dgs{idgs};
    
    
end