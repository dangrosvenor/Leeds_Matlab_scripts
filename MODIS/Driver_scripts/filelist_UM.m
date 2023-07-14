clear fileUM2

idat=1;

for i=1:99
    flag{i}='';
end

fileUM2{idat} = 'xkqkh_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkh) 100cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM2{idat} = 'xkqkj_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkj) 400cm^{-3}';  pole_lat=70; pole_lon=278;idat=idat+1;
fileUM2{idat} = 'xkqkk_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkk) 400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM2{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM2{idat} = 'xkqko_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqko) 100cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM2{idat} = 'xkqkq_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkq) 100cm^{-3} No cloud-scheme'; pole_lat=70; pole_lon=278;idat=idat+1;
fileUM2{idat} = 'xkqkr_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkr) 1000cm^{-3} No cloud-scheme';pole_lat=70; pole_lon=278; idat=idat+1;
fileUM2{idat} = 'xkqkf_qL_qR_.pp.nc.mat'; labs_UM(idat).l = '(xkqkf) Old-mphys'; flag{idat}='load_mat'; fileUM_rho{idat} = 'xkqkf_rho_.pp.nc'; pole_lat=70; pole_lon=278; idat=idat+1;


for idat=1:length(fileUM2)
    
    
end
