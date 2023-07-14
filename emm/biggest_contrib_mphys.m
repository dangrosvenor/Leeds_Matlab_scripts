%18/09/07 script to find the biggest contributer from emm mphys for a particular height

clear totmass totnum

idir=1;

h=10.25;
h=10.4;

ih=findheight(Zmphys/1000,h);
fprintf(1,'\n****Height used is %f *****\n',Zmphys(ih)/1000);
for im=3:size(mphys(idir).mass,3)               %%%%%%%%% first two columns are time and height so ignore
    totmass(:,im-2)=sum(mphys(idir).mass(:,:,im),2);
    totnum(:,im-2)=sum(mphys(idir).nums(:,:,im),2);
end



    
% [maxval imax]=max(totmass);
 ntim=0.5/(0.2/60);  %no. timesteps - divided by this in early version
% totmass(imax)*1000*ntim %mulitply by this as there was no need to divide by it in the first place 
%sum of this should give total mass (g/kg) or num created

%dgs={'ihomog','ievap_num','icond','ievap_size','iauto','ievap_zbase','ievap_ctop','iauto_res','iauto_purge'};60)

figure
plot(1000*ntim*totmass(:,12),Zmphys/1000);