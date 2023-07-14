function Zs=wrf_radar_Zs_1(N0_exp,lam_exp,am_s,rho_w,bm_s)
%bid to save memory by only temporarily creating arrays in these functions
%function Zs=wrf_radar_Zs_1(N0_exp,lam_exp,am_s,rho_w,bm_s)

Zs=N0_exp.*(6.*am_s./pi./rho_w).^2./lam_exp.^(1+0+2.*bm_s) ...
    .*gamma(1+0+2.*bm_s);