function [Ty2,RHfp,RHhumi,RHfp_i,RHhumi_i,LWC_fp,Tad_fp,LWC_humi...
    ,Tad_humi,qv_wave_fp,qv_wave_humi,py]...
    =Twave(y2,h,pot,P,Zsound,qfp_sound,qhumi_sound,qvfp_mean,qvhumi_mean,iuse_mean_extrap...
    ,iuse_constant_RH_above_sounding,RH,svp_method_liq,svp_method_ice)


%iuse_mean_extrap=1; %whether to extrapolate qv to alitudes we don't have data for by assuming a
            %constant value

f=1e6*28.97/18;

potx0=interp1(h,pot,y2(1)); %starting potemp
%now work out temp at each displacement in y
py = interp1(h,P,y2); %pressure for all displacements
Ty2 = potx0 ./ (1000e2./py).^0.286 - 273.15; %assuming constant potemp

qvfp_x0=interp1(h,qfp_sound,y2(1)); %starting qv
if iuse_mean_extrap==1 & y2(1)>Zsound(end)
    qvfp_x0=qvfp_mean;
end

qvhumi_x0=interp1(h,qhumi_sound,y2(1)); %starting qv
if iuse_mean_extrap==1 & y2(1)>Zsound(end)
    qvhumi_x0=qvhumi_mean;
end

if iuse_constant_RH_above_sounding==1 & y2(1)>Zsound(end)
    qsat_x0 = satvappress(273.15+Ty2(1),svp_method_liq,'liq',py(1),1)/f;
    qvfp_x0 = RH*qsat_x0;
    qvhumi_x0 = RH*qsat_x0;
end

qsatMR = (SatVapPress(Ty2+273.15,svp_method_liq,'liq',py,1)/f); 
RHfp = 100 * qvfp_x0 ./ qsatMR;
RHhumi = 100 * qvhumi_x0 ./ qsatMR;
qv_wave_fp = qvfp_x0*ones([1 length(qsatMR)]);
qv_wave_humi = qvhumi_x0*ones([1 length(qsatMR)]);

RHfp_i = 100 * qvfp_x0 ./ (SatVapPress(Ty2+273.15,svp_method_ice,'ice',py,1)/f);
RHhumi_i = 100 * qvhumi_x0 ./ (SatVapPress(Ty2+273.15,svp_method_ice,'ice',py,1)/f);

icloud_base = find(qsatMR<=qvfp_x0); %find all the saturated positions
icloud_base2 = find(qsatMR<=qvhumi_x0); %find all the saturated positions

if length(icloud_base)>0
    cbP_fp = py(icloud_base(1)); %cloud base
    cbT_fp = Ty2(icloud_base(1))+273.15;    
end
LWC_fp=NaN*ones(1,length(Ty2)); %set the vector to NaN to start with
%Tad_fp=NaN*ones(1,length(Ty2)); %set to NaN to start with
Tad_fp=Ty2; %set to NaN to start with

if length(icloud_base2)>0
    cbP_humi = py(icloud_base2(1)); %cloud base
    cbT_humi = Ty2(icloud_base2(1))+273.15;
end
LWC_humi=NaN*ones(1,length(Ty2)); %set the vector to zero to start with
%Tad_humi=NaN*ones(1,length(Ty2)); %set to NaN to start with
Tad_humi=Ty2;

for i=1:length(icloud_base)
    [LWC_fp(icloud_base(i)),Tad_fp(icloud_base(i))]...
        =adLWC_PaulLawson_simple(cbP_fp,cbT_fp,py(icloud_base(i)));
    
   rho = density(py(icloud_base(i)),273.15+Tad_humi(icloud_base(i))); 
   qv_wave_fp(icloud_base(i)) = qv_wave_fp(1) - LWC_fp(icloud_base(i))/1000/rho;
    
end

for i=1:length(icloud_base2)
    [LWC_humi(icloud_base2(i)),Tad_humi(icloud_base2(i))]...
        =adLWC_PaulLawson_simple(cbP_humi,cbT_humi,py(icloud_base2(i)));
    
%     qv_wave_humi(icloud_base2(i)) = SatVapPress(Tad_humi(icloud_base2(i))+273.15,'goff','ice',py(icloud_base2(i)),1)/f;
      rho = density(py(icloud_base2(i)),273.15+Tad_humi(icloud_base2(i)));
      qv_wave_humi(icloud_base2(i)) = qv_wave_humi(1) - LWC_humi(icloud_base2(i))/1000/rho;
end



