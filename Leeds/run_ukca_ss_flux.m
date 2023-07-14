u=[1:50];
u=[4.1 5 10];
%for iu=1:length(u)
%[aer_mas_primss,aer_num_primss,aer_numflux_primss,aer_massflux_primss,drlo,drmid,drhi] = ukca_prim_ss(u(iu));
[numflux_accum, numflux_coarse, massflux_accum, massflux_coarse,drlo,drmid,drhi] = ukca_prim_ss(u)
%numflux_accum(iu)=aer_numflux_primss(:,:,:,3);
%numflux_coarse(iu)=aer_numflux_primss(:,:,:,4);
%massflux_accum(iu)=aer_massflux_primss(:,:,:,3);
%massflux_coarse(iu)=aer_massflux_primss(:,:,:,4);
%end

numflux_accum
numflux_coarse
massflux_accum
massflux_coarse

figure
plot(u,numflux_accum); grid on; hold on
plot(u,numflux_coarse,'k')
xlabel('U10 (m s^{-1})');
ylabel('Total SS number flux from surface (# m^{-2} s{-1})');
set(gca,'yscale','log');
legend({'Accumulation mode','Coarse mode'});

figure
plot(u,massflux_accum); grid on; hold on
plot(u,massflux_coarse,'k')
xlabel('U10 (m s^{-1})');
ylabel('Total SS mass flux from surface (kg m^{-2} s{-1})');
set(gca,'yscale','log');
legend({'Accumulation mode','Coarse mode'});

drlo(1)*1e6
drhi(end)*1e6
drmid(1)*1e6
drmid(end)*1e6

