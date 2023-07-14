function heff=Houghton_find_heff_from_streamline_del(del,U,gd,H0_Houghton,HM,side)

range=[1 HM];
   

heff = fzero(@Houghton_find_heff_from_streamline_del_fun,range,[],del,U,gd,H0_Houghton,side);

