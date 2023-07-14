function y=Houghton_find_heff_from_streamline_del_fun(h,del,U,gd,H0_Houghton,side)



thi=mountain_Houghton_solve_thi_for_given_h(h,U,H0_Houghton,gd,side);
y = H0_Houghton - h + del - thi; %should be zero since thi=H0_Hougton-h+del





