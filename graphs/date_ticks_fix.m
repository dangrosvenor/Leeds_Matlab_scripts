[temp,idat_min] = min(start_time2);
[temp,idat_max] = max(end_time);

       [Y,M,D,H,MM] = datevec(xdat(idat_min).x);
       [H1,imin] = min(H);
       H=H1:H1+6;
       istart_temp=find(rem(H,6)==0);
       start_val=H(istart_temp(1));
       Htime = datenum(Y(imin),M(imin),D(imin),start_val,0,0); 
            
              
       xtickvals_set = [Htime : 6/24 : max(xdat(idat_max).x)];
        
        [Y,M,D,H,MM] = datevec(xtickvals_set);
        month = datestr(xtickvals_set,'mmm');
        time_str144 = datestr(xtickvals_set,'HH');
        
        
        
        for ixtick_lab=1:length(xtickvals_set)
            xticklabs_set{ixtick_lab} = [time_str144(ixtick_lab,:)];
        end
        
       istart_temp=find(H==0);
        for ixtick_lab=[istart_temp:4:length(xtickvals_set)]
            xticklabs_set{ixtick_lab} = [num2str(D(ixtick_lab)) '-' month(ixtick_lab,:)];
        end 
        
     iset_xticks = 1;
     iset_xticklabs=1;