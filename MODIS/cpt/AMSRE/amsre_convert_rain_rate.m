function [tlwp_out] = amsre_convert_rain_rate(lwp, rain, sst )
%Calculates the total liquid water path (LWP+RWP) from amsre LWP and rain
%rate.
%rain is the rain rate
%This is likely to be an overestimate (I think) and therefore should be an
%upper limit of the TLWP.
%This is because if the amount of rain is not diagnosed then there would be too much attentuation, which would get interpreted as 
%extra condensate (since rain distributions attenuate twice as efficiently
%as cloud DSDs according to Lebsock, JGR, 2014 (Section 3.4.4).


nx=size(lwp,1);
ny=size(lwp,2);
nt=size(lwp,3);
nasc=size(lwp,4);

tlwp_out = NaN*ones(size(lwp));

max_t_chunk = 200; %split into chuncks of 200 in the time dimension
tinds = [1:max_t_chunk:nt];

it_start=1;
it_end=tinds(2)-1;
for it = 1: length(tinds)-1
    [tlwp_out(:,:,it_start:it_end,:)] = amsre_convert_rain_rate_partial(lwp(:,:,it_start:it_end,:), rain(:,:,it_start:it_end,:), sst(:,:,it_start:it_end,:) );
    it_start = it_start+max_t_chunk;
    it_end = min(it_end+max_t_chunk, nt);
    
end



function [tlwp_out] = amsre_convert_rain_rate_partial(lwp, rain, sst_in )


conv=1/3600.0;
rhow=1000.0;
Rair=287.05;
Cp=1005.0;
pref=1e+05;
kappa=Rair/Cp;


%Abel-Boutle DSD;
%;
Na=0.22;
Nb=2.2;
aa=0.0;


%m-D;
%;
a=pi*rhow/6.0;
b=3.0;

%vt-D;
%;
c=362.0;
d=0.65;

alpha=0.18;
const=0.16;

h0=0.46;
h1=5.16;

%satelite data;
%;

tlwp_out = NaN*ones(size(lwp));
hh=NaN*ones(size(lwp));
rwp1 = NaN*ones(size(lwp));
%hh=NaN*ones([nx ny nt]);
%rwp1 = NaN*ones([nx ny nt]);

%mlwp=fltarr(np,ny);

%lons=0.25*(findgen(nx)+1)- 0.125;
%lats=0.25*(findgen(ny)+1)-90.125;


%Duplicate the sst for both the day and night overpasses
sst = repmat(sst_in,[1 1 1 2]);



%for i=0,nfile-1;

%     file_name=amsr_data(i);
% 
%     read_amsr_day_v7, file_name, time1, sst1, wind_lf1, wind_mf1, vapor1, cloud1, rain1;
% 
%     file_name=amsr_3d3(i);
%     read_amsr_averaged_v7, file_name, sst_3d3, wind_lf_3d3, wind_mf_3d3, vapor_3d3, cloud_3d3, rain_3d3;
% 
%     tmp1=sst_3d3;
    %% tmp1=sst1;

%    tmp1 is the sst field from the 3d3 files (3-day means?)
%    rain are the rain rates

% ----------------------------------------------------------------------
% N.B. - Kalli used a threshold of 200 here to signify bad data - in my
% version of the read file I have set all of these to NaN already.
% Plus the arrays are initialized as NaN to start with.


%    rain=rain1(*,*,0);
    %ix=where(lwp lt 200 and rain lt 200,complement=ixc);
%    ix=find(sst < 200 & rain < 200);

    ix=~isnan(sst) & ~isnan(rain);

%    ixc = [1:prod(size(lwp))]; ixc(ix) = []; %the complement of ix
%    ixc = setdiff([1:prod(size(lwp))],ix)

    %AMSR assumes rain is distributed over a;
    % a height (hh) which is a function of the SST;


    hh(ix)=h0+const*sst(ix);            %[km];
    
    
%     i0 = where(tmp1 lt 0.0 and tmp2 lt 200)
%     i1 = where(tmp1 gt 30.0 and tmp1 lt 200 and tmp2 lt 200)
%     i2 = where(tmp1 ge 200 or tmp2 ge 200)

    i0 = find(sst < 0.0 & isnan(rain)==0);
    i1 = find(sst > 30.0 & isnan(sst)==0 & isnan(rain)==0);
    i2 = find(isnan(sst)==1 | isnan(rain)==1);

    hh(i0)=h0; %in km
    hh(i1)=h1;
    hh(i2)=NaN;

    %use the DSD to convert the AMSR rain rate to a rwc;
    %;
%    rain_flux =  conv .* rain;

    rwc =  a.*Na .* gamma(b+aa+1) .* (conv.*rain./(a.*c.*Na.*gamma(b+d+aa+1))).^((b+aa+1-Nb)./(b+d+aa+1-Nb));
    %rwc =  rain_flux/7.0;

    %assume the rain water is uniform over height hh;
    % to get a rwp;

%    rwp1 = fltarr(nx,ny);

    rwp1(ix) = 1e+03.*hh(ix).*rwc(ix);          % [kg/m^2];
%    rwp1(ixc)=-99.0;

%    lwp1=fltarr(nx,ny);
%    lwp1 = NaN*ones(size(lwp));
%    lwp = cloud1(*,*,0);


    %set -ve to zero;
    %;
    %% ix=where(lwp lt 0.0);
    %% lwp(ix)=0.0;

    %rain + cloud is the total liquid water path;
    %;
%    ix=where(sst lt 200 and rain lt 200 and sst lt 200 and sst ge -0.05,complement=ixc);
    ix=find(isnan(sst)==0 & isnan(rain)==0 & isnan(lwp)==0 & lwp >= -0.05);    
%    ixc = [1:prod(size(sst))]; ixc(ix) = []; %the complement of ix
    
    tlwp_out(ix) = rwp1(ix) + lwp(ix);
%    tlwp_out(ixc)=-99.0;

%    lwp(2*i,*,*)=lwp1;

%calculates the mean vs lat
%     for j=0, ny-1;
%         tmp=lwp1(*,j);
%         ix=where(tmp gt -99.0,countx);
%         if countx ne 0;
%             mlwp(2*i,j)=mean(tmp(ix));
%         end;
%     end;
'';

