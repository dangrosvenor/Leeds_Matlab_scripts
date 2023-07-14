file_um = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-bj045/umnsaa_pc030.pp.nc';
nc_fix = netcdf(file_um);
file_um_bug = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-bj048/umnsaa_pc030.pp.nc';
nc_bug = netcdf(file_um_bug);

iy=150;

x=nc_bug{1}(:);
y=nc_bug{2}(:);
z=nc_bug{3}(:)/1e3; %km

%order is 't, hybrid_ht, y, x'
nd_bug_ALL=nc_bug{10}(:,:,:,:)/1e6; %per cc
nd_fix_ALL=nc_fix{10}(:,:,:,:)/1e6;
cf_bug_ALL=nc_bug{9}(:,:,:,:);
cf_fix_ALL=nc_fix{9}(:,:,:,:);

inan=find(nd_bug_ALL<1); nd_bug_ALL(inan)=NaN;
inan=find(nd_fix_ALL<1); nd_fix_ALL(inan)=NaN;
nd_bug_ALL = nd_bug_ALL ./ cf_bug_ALL; %These are all either zero or 1 anyway...
nd_fix_ALL = nd_fix_ALL ./ cf_fix_ALL;

mean_bug = meanNoNan(nd_bug_ALL(:),1);
mean_fix = meanNoNan(nd_fix_ALL(:),1);



nd_bug=nd_bug_ALL(1,:,iy,:); %per cc
nd_fix=nd_fix_ALL(1,:,iy,:);

figure
pcolor(x,z,squeeze(nd_bug)); shading flat; colorbar
set(gca,'ylim',[0 10]);
set(gca,'xlim',[179 183]);
caxis([0 300]);
ylabel('Height (km)');
xlabel('Longitude (degrees)');
title('N_d (cm^{-3}) WITH bug');


figure
pcolor(x,z,squeeze(nd_fix)); shading flat; colorbar
set(gca,'ylim',[0 10]);
set(gca,'xlim',[179 183]);
caxis([0 300]);
ylabel('Height (km)');
xlabel('Longitude (degrees)');
title('N_d (cm^{-3}) WITHOUT bug');




iz=find(z>10);
iz=iz(1)-1;
nd_bug_10km = nd_bug_ALL(:,iz,:,:); %per cc
nd_fix_10km = nd_fix_ALL(:,iz,:,:);
mean_10km_bug = meanNoNan(nd_bug_10km(:),1)
mean_10km_fix = meanNoNan(nd_fix_10km(:),1)

%mean profile
prof_10km_bug = meanNoNan(meanNoNan(meanNoNan(nd_bug_ALL,1),3),2);
prof_10km_fix = meanNoNan(meanNoNan(meanNoNan(nd_fix_ALL,1),3),2);
figure
plot(prof_10km_bug,z,'r','linewidth',3);
hold on
plot(prof_10km_fix,z,'b--','linewidth',3);
set(gca,'ylim',[0 18]);
leg_str={'Bug','fix'};
legend(leg_str);
ylabel('Height (km)');
xlabel('N_d (cm^{-3})');

prc_diff = 100*( (prof_10km_bug - prof_10km_fix) ./ prof_10km_fix );
figure
plot(prc_diff,z,'r','linewidth',3);
set(gca,'ylim',[0 18]);
ylabel('Height (km)');
xlabel('Percentage difference (bug vs fix; %)');

abs_diff = prof_10km_bug - prof_10km_fix;
figure
plot(abs_diff,z,'r','linewidth',3);
set(gca,'ylim',[0 18]);
ylabel('Height (km)');
xlabel('Difference (bug vs fix; cm^{-3})');



