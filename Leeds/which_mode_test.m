f=[1 0.99 0.5 0.1 0.01 0.001];
methods={'UM pre bug-fix','UM','Robin','Robin2','Dan Robin','Dan Robin2','lognormal dN/lnD ratio','midway all or nothing'}
%methods={'UM'}

for im=1:length(methods)
method = methods{im};

for i=1:length(f)
    
  ff= f(i);    
  f_init = 0.5;
  [dm1(im,i), dm2(im,i), dn1(im,i), dn2(im,i), dm_ratio_accum(im,i), dn_ratio_accum(im,i),r1_new(im,i),r2_new(im,i),rm(im,i),r_ratio(im,i),r1(im,i),r2(im,i)] = which_mode(f_init*4.5e-9, f_init*ff*3.8e8,method);
  %N.B. - multiplying the mass and number by 0.5 to simulate effect of
  %starting with half the namelist mass and number and adding aerosol


end


end