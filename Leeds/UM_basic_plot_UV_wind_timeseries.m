%plot sqrt(U.^2+V.^2) timseries (domain mean)


file_UM = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-al522/u-al522_U_wind_10m_VOCALS_1p0_L70_ukv_.pp.nc.mat';
load(file_UM,'U_wind_10m');
u10_timser = meanNoNan(meanNoNan(U_wind_10m,2),2);

file_UM = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-al522/u-al522_V_wind_10m_VOCALS_1p0_L70_ukv_.pp.nc.mat';
load(file_UM,'V_wind_10m');
v10_timser = meanNoNan(meanNoNan(V_wind_10m,2),2);

WS_timser = sqrt(u10_timser.^2+v10_timser.^2);

figure; clear leg_cell; ileg=1;
% plot(WS_timser,'bo-'); hold on; 
% leg_cell{ileg}='Using domain mean U and V'; ileg=ileg+1;

WS = sqrt( U_wind_10m.^2 + V_wind_10m(:,1:end-1,:).^2 ); %Ignoring the staggering here for now
WS_timser2 = meanNoNan(meanNoNan(WS,2),2);

plot(WS_timser2,'rx-'); hold on
leg_cell{ileg}='Mean WS'; ileg=ileg+1;

%From Paul - I guess the other thing is to average the winds according to the form of the flux-wind relation.
%So it looks exponential ish. If you average the exp(wind) and then log it what wind speed does that give?

eWS=exp(WS);
WS_timser_exp = log ( meanNoNan(meanNoNan(eWS,2),2) );
plot(WS_timser_exp,'g^-'); hold on
leg_cell{ileg}='ln(mean(exp(WS))'; ileg=ileg+1;

legend(leg_cell);
xlabel('Time index');
ylabel('Wind speed (m s^{-1})');

%Mean ignoring first few indices
mean01=meanNoNan(WS_timser(3:end),1);
mean02=meanNoNan(WS_timser2(3:end),1);
fprintf(1,'\nMean of WS_timser=%f; Mean of WS_timser2=%f\n',mean01,mean02);

std01=std(WS_timser(3:end),1);
std02=std(WS_timser2(3:end),1);
fprintf(1,'\nStd dev of WS_timser=%f; Std dev of WS_timser2=%f\n',std01,std02);

it=8; %time of max WS after initial peak
WS_01 = squeeze(WS(8,:,:));
[numflux_accum, numflux_coarse, massflux_accum, massflux_coarse,drlo,drmid,drhi] = ukca_prim_ss(WS_01);
numflux_accum_nonlinear = meanNoNan(numflux_accum(:),1);
[numflux_accum_mean, numflux_coarse, massflux_accum, massflux_coarse,drlo,drmid,drhi] = ukca_prim_ss(WS_timser2(it));

%If have a flux of 2e5 particles per m2 per s and a 1000m BL height then
%would replenish
%For U=6 m/s looks like would expect a flux of around F=2e5 from Martensson
%(2003) Fig. 10. Since we have windspeeds that are more like 4 m/s on
%average, this is perhaps an overestimate.
F=2e5;
dN_dhr = F*1e-6*3600/1000; %number per cc per hour

%With F=2e5 gives only 0.7 per cc per hour. So, quite slow...





