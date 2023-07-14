function [n_time,n_ice,n_drops,time,r_ice,c_ice,r_drops,c_drops,...
        RH,c_ice_tot,sum_c,c_end,sum_cy1,sum_cy2,arvot]=read_data_ice(dirpath,out_folder)
%WRITE statement in cmodel is as follows:
% WRITE(IOUT,459)time,(D_Y1(I),I=1,NFBins+1),(C_Y1(i),i=1,NFBins), &
%           &     (r_Y2(I),I=1,NY2+1),(C_Y2(i),i=1,NY2),RHI,C_ice_tot,sum(c(NABins+1:2*NABins)),C(InC(IH2Ogi,NABins)), &
%           &  sum(c_Y1(NY1+1:int(2)*NY1)),sum(c_Y2(NY2+1:int(2)*NY2))

          
          
NFbins = 60; %no. ice bins
NY2=60; %no. bins for moving centre droplets
Nfields = NFbins*2+1 + NY2*2+1 + 7;

max_size=1e3;

filename=[dirpath out_folder 'iceout001.out'];

FID=fopen(filename,'rt');
[arvot,COUNT] = fscanf(FID,'%f',[Nfields,Inf]);

%ice=dlmread([dirpath out_folder 'iceout001.out'],' ',[0 0 max_size-1 Nfields-1]);
%ice=dlmread([dirpath out_folder 'iceout001.out'],' ',[0 0 Nfields-1 max_size-1]);
arvot=arvot';
 
  time=arvot(:,1);                        n=1;
  r_ice=arvot(:,n+1:n+NFbins+1);          n=n+NFbins+1;
  c_ice=arvot(:,n+1:n+NFbins);            n=n+NFbins;
  r_drops=arvot(:,n+1:n+NY2+1);           n=n+NY2+1;
  c_drops=arvot(:,n+1:n+NY2);             n=n+NY2;
  RH=arvot(:,n+1);                        n=n+1;         
  c_ice_tot=arvot(:,n+1);                 n=n+1;
  sum_c=arvot(:,n+1);                     n=n+1;
  c_end=arvot(:,n+1);                     n=n+1;
  sum_cy1=arvot(:,n+1);                   n=n+1;
  sum_cy2=arvot(:,n+1);                   n=n+1;
  
  n_time=length(time);
  n_ice =size(c_ice,1);
  n_drops =size(c_drops,1);

%--------------------------------------------------------------------------