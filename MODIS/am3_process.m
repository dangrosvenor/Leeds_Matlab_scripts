% plot_global_maps  (use to shorcut to it)
icf_low_only=0; %flag to say whether we want to screen the CF for LWP using only low clouds (pressure
%screening), or for all clouds
i_ice_screen=1; %flag to screen for the presence of any ice - if a layer has ice then the cloud 
%fraction in that layer is ignored and not used for determining whether we have a CF>0.8
%in a column, which is used for determining whether to include the LWP value for that column

min_cf=0.01;

low_liq_thresh = 0.05e-3; %kg/kg
low_cf_thresh = [0.01 1]; %will be excluding isolated clouds for low CF values
low_cf_thresh = [0.8 1.01]; %
%low_cf_thresh = [0.0 0.8]; %
thresh_Nlevs = 0; %minimum number of vertical model levels required for a vertically
thresh_P = [750 1080]*1e2;
thresh_Nd = 5; %per cc

%have divided by cloud fractions to get in-cloud averages - remove values when
%have low CF, as will be inflated by zero cloud-fractions. So we insure a
%minimum CF here (=min_cf)
ii=find(low_cf_thresh<min_cf);
low_cf_thresh(ii)=min_cf;

%indices for points that DO NOT have the desired cloud fraction
ilow_cf = find(~ (am3_cf>low_cf_thresh(1) & am3_cf<=low_cf_thresh(2)) );

%points to remove that are out of our pressure range
ip = find( ~ ( am3_pfull>=thresh_P(1) & am3_pfull<thresh_P(2) )  );

iNd = find( 1e-6*am3_drop < thresh_Nd ); %to remove for a few low Nd regions with reasonable LWC
%e.g. up to 0.05 g/kg, but with only 2 per cc droplets - likely rain I
%guess, or evaporating old cloud. Likely to have a low tau anyway, so MODIS
%likely wouldn't detect them

%combination of the two sets of indices to be removed with no repetitions
%(can use union or unique)
iremove = union(ilow_cf,[ilow_cf; iNd]);
%iremove = unique([ilow_cf; ip; iNd]);

am3_drop2 = 1e-6 * am3_drop; %cm3 - already done the /kg to /m3 conversion
am3_liq2 = am3_liq;
am3_pfull2 = am3_pfull;

am3_CF_screened = am3_cf;
am3_cf2 = am3_cf;


if icf_low_only==1
    %for CF want to screen for LWC, Nd and Pressure, but not CF
    am3_CF_screened(ip)=NaN; %set cloud fraction to NaN when outside of required pressure range
end


if i_ice_screen==1
    i_ice = find(am3_iwc_av>0.0001e-3);
    am3_CF_screened(i_ice)=NaN; %set cloud fraction to NaN when outside of required range
end

i_remove_CF = find(am3_liq<low_liq_thresh);
%am3_CF_screened(i_remove_CF)=0; %set cloud fraction to zero when there is little LWC

%when have v.low cf will assume LWC=0 and ignore CF (or set to zero)
i_remove_CF = find(am3_cf<min_cf); %remove CF when CF is near zero as are dividing by 0
am3_CF_screened(i_remove_CF)=0; %set cloud fraction to zero when there is little CF - perhaps not necessary
%since CF is already low!! - Although want to remove even low CFs when the
%LWC is low - am assuming that don't know LWC when CF<0.01 - are replacing
%e.g. CF=0.009 with CF=0 since the cloud is likely undetectable
%do we need to consider the effect of MODIS pixel size?

i_remove_CF = find(am3_drop2<thresh_Nd);
%am3_CF_screened(i_remove_CF)=0; %set cloud fraction to zero when there are few droplets

am3_CF_max_screened = squeeze(max(am3_CF_screened,[],2));
am3_CF_mean_screened = squeeze(meanNoNan(am3_CF_screened,2));



%here we screen for pressure, Nd and cloud fraction (at each height)
am3_drop2(iremove)=NaN;
am3_liq2(iremove)=NaN;
am3_cf2(iremove)=NaN;
am3_pfull2(iremove)=NaN;

%am3_CF_min is therefore the min CF of the layers that meet our above
%criteria - doesn't matter that we have removed some CFs 
%for the purpose of including only layers where
%the min CF is above low_cf_thresh(1), since the ones we removed were all below
%low_cf_thresh(1)
%am3_CF_min = squeeze(min(am3_cf2,[],2));
am3_CF_min = squeeze(max(am3_CF_screened,[],2));
%now we have removed the height dimension - use to screen for lwp
%actually we know the min cannot be below 0.8! So could just check for columns
%where have some numbers other than NaN - this is akin to doing that
iadd_lwp = find(am3_CF_min>=low_cf_thresh(1)); %use LWP only when we have at least thresh CF somewhere
%in the column
%start with an array of NaNs
am3_lwp_minCF = NaN*ones(size(am3_lwp_mean));
%add the desired ones
am3_lwp_minCF(iadd_lwp)=am3_lwp_mean(iadd_lwp);
%this will now be the cell average LWP, rather than the cloud-average
%is probably best to compare this to the MODIS cell average value
am3_lwc_minCF = NaN*ones(size(am3_liq));
%add the desired ones
[T,I,J] = ind2sub(size(am3_lwc_minCF),iadd_lwp);
%K=repmat([1:size(am3_liq,2)],[size(am3_liq,1) 1 size(am3_liq,3) size(am3_liq,4)]);
K=repmat([1:size(am3_liq,2)],[length(I) 1]);
K=K';
T2=repmat(T,[1 size(am3_liq,2)]); T2=T2';
I2=repmat(I,[1 size(am3_liq,2)]); I2=I2';
J2=repmat(J,[1 size(am3_liq,2)]); J2=J2';
lwc_inds = sub2ind(size(am3_liq),T2(:),K(:),I2(:),J2(:));

am3_lwc_minCF(lwc_inds)=am3_liq_av(lwc_inds);
am3_lwc_minCF2=am3_lwc_minCF;
%remove lwc if we are out of the CF criteria
am3_lwc_minCF2(ilow_cf)=NaN;

i_remove_lwc = find(am3_lwc_minCF2<low_liq_thresh);
%for cloud depth probably just want to use layers when CF>0.8
am3_pfull_ice_screen = am3_pfull;
am3_pfull_ice_screen(ilow_cf) = NaN;
am3_pfull_ice_screen(i_ice) = NaN; %remove ice points
am3_pfull_ice_screen(iNd) = NaN;
am3_pfull_ice_screen(i_remove_lwc) = NaN;
%am3_pfull_ice_screen(ip) = NaN;

am3_Nd_ice_screen = 1e-6*am3_drop;
am3_Nd_ice_screen(ilow_cf) = NaN;
am3_Nd_ice_screen(i_ice) = NaN; %remove ice points
am3_Nd_ice_screen(iNd) = NaN;
am3_Nd_ice_screen(i_remove_lwc) = NaN;
am3_Nd_max_screen = squeeze( max(am3_Nd_ice_screen,[],2) );

%am3_pfull_ice_screen(ip) = NaN;




%end of new ice-cloud fraction screening stuff



i_low_liq2 = find(am3_liq2<low_liq_thresh);
am3_drop2(i_low_liq2)=NaN;
am3_pfull2(i_low_liq2)=NaN;
am3_liq2(i_low_liq2)=NaN;

%am3_drop_ice_screen = NaN*ones(size(am3_drop));
%am3_drop_ice_screen(lwc_inds) = 1e-6*am3_drop(lwc_inds);
%am3_drop_ice_screen(ilow_cf) = NaN;  % 

am3_drop_ice_screen = 1e-6*am3_drop;
am3_drop_ice_screen(ilow_cf) = NaN;  %need to do this because of the divide by zero CF problem 

am3_Nd_maxALL = squeeze(max(am3_drop2,[],2));
am3_Nd_maxALL_ice_screen = squeeze(max(am3_drop_ice_screen,[],2));





% am3_cf_thresh = 0.8;
% icf = find( ~(am3_cf>=am3_cf_thresh & am3_liq2./am3_cf>=am3_liq_thresh) );
% %create an array of the model levels (1=TOA, 48=ground), so that we
% %can remove the points that don't meet the cloud criteria. Then the lowest
% %number for the column of each cell will denote the highest model level
% %that meets our cloud criteria
% s_am3=size(am3_cf);
% am3_icf = repmat([1:s_am3]',[1 s_am3(2) s_am3(3)]);
% am3_icf_orig = am3_icf;
% am3_icf(icf)=s_am3(1); %set to min - the surface
% am3_maxz = min(am3_icf,[],1); %the surface is actually at iz=48 so we want the min index
% 
% maxz_lin = am3_maxz(:);
% am3_maxzp = am3_pfull(maxz_lin);
% am3_maxzp = reshape(am3_maxzp,size(am3_maxz));

%could take the Nd at the point with the highest LWC in each profile
%or could take the mean Nd over the profile - Nd seems to vary quite a lot
%in the model profile I looked at - need to look at more

if i_ice_screen==1
    am3_pfull2 = am3_pfull_ice_screen;
end
   

am3_liq2 = permute(am3_liq2,[2 1 3 4]);
am3_drop2 = permute(am3_drop2,[2 1 3 4]);
am3_drop_ice_screen = permute(am3_drop_ice_screen,[2 1 3 4]);
am3_cf2 = permute(am3_cf,[2 1 3 4]);
h_half2 = permute(h_half,[2 1 3 4]);
am3_pfull2 = permute(am3_pfull2,[2 1 3 4]);



s_am3_low = size(am3_cf2);

[max_liq,imax_liq] = max(am3_liq2,[],1);
max_liq = squeeze(max_liq);

%do the same for the cloud top and base to calculate cloud depth
[am3_CTP,iCTP] = min(am3_pfull2,[],1); %find the minimum pressure of present LWC (cloud top)
%[am3_CBP,iCBP] = max(am3_pfull2,[],1); %find the maximum pressure of present LWC (cloud base)
%have to be concerned about multiple cloud layers, though
am3_CTP = squeeze(am3_CTP);
%am3_CBP = squeeze(am3_CBP);

%use the highest continuous cloud layer - so cloud top is the min pressure
%with LWC as above. Have put NaNs at locations where there is no cloud
%meeting our criteria - so look for continuous layers of non-NaN values.
inan=isnan(am3_pfull2);
inan(end+1,:,:,:)=1; %add an a column to the height dimension to make the diff array the same size
%as the original - this way if the height index closest to the surface is not NaN (inan=0) then diff
%will equal 1 for the surface index and will report the lowest level as
%the cloud base height -  diff does X(2:end)-X(1:end-1) will

d=diff(inan,1,1);
ii=find(d==1);

iCBP2 = NaN*ones(size(d));
[K,T,I,J]=ind2sub(size(d),ii);
iCBP2(ii) = K;

iCBP = min(iCBP2,[],1); %find the min value since this corresponds to the highest altitude level
iremove2 = find(isnan(iCBP)==1); %points where have no cloud
iCBP(iremove2)=1; %cannot have an index as NaN or will cause an error
%will remove these values later



%imax_liq now gives a matrix containing the z-index of the max LWC for each
%position and time

%convert these into a linear index for the the whole of the big matrix
% use IND = SUB2IND(SIZ,I,J)
% will use imax_liq(:) - i.e. a vector containing a list of all the
% z-indices. SIZ will be [NZ NT*NI*NJ]. I will be the NZ indices. J will be
% all of the other indices. i.e. 1:NT*NI*NJ
imax_lin = sub2ind([s_am3_low(1) s_am3_low(2)*s_am3_low(3)*s_am3_low(4)],imax_liq(:),[1:s_am3_low(2)*s_am3_low(3)*s_am3_low(4)]');
%max_test = reshape(am3_liq2(imax_lin),s_am3_low(2),s_am3_low(3),s_am3_low(4)); %are the same

%for heights we will use the h_half2 array which is one level larger than
%h_full in the vertical. Level N = 0m, level N-1 = the top of the first layer
%in h_half2 (h_half2(N) is between the h_full(N) and h_full(N-1) )
%so will find the linear indices for the larger matrix h_half2 and will use
%the height indices +1 for CTH & the index as is for CBH
imax_lin_CTH = sub2ind([s_am3_low(1)+1 s_am3_low(2)*s_am3_low(3)*s_am3_low(4)],iCTP(:),[1:s_am3_low(2)*s_am3_low(3)*s_am3_low(4)]');
imax_lin_CBH = sub2ind([s_am3_low(1)+1 s_am3_low(2)*s_am3_low(3)*s_am3_low(4)],iCBP(:)+1,[1:s_am3_low(2)*s_am3_low(3)*s_am3_low(4)]');
imax_lin_CBP = sub2ind([s_am3_low(1) s_am3_low(2)*s_am3_low(3)*s_am3_low(4)],iCBP(:),[1:s_am3_low(2)*s_am3_low(3)*s_am3_low(4)]');

am3_Nd_maxliq = reshape(am3_drop2(imax_lin),s_am3_low(2),s_am3_low(3),s_am3_low(4)); %
%have already converted to per cc
am3_CF_maxliq = reshape(am3_cf2(imax_lin),s_am3_low(2),s_am3_low(3),s_am3_low(4)); 
am3_liq_maxliq = reshape(am3_liq2(imax_lin),s_am3_low(2),s_am3_low(3),s_am3_low(4)); 

am3_CTH = reshape(h_half2(imax_lin_CTH),s_am3_low(2),s_am3_low(3),s_am3_low(4)); %have now lost the height dimension
am3_CBH = reshape(h_half2(imax_lin_CBH),s_am3_low(2),s_am3_low(3),s_am3_low(4));
am3_CBP = reshape(am3_pfull2(imax_lin_CBP),s_am3_low(2),s_am3_low(3),s_am3_low(4));

%script to take all the iCBP and iCTP indices and make a 3D LWC array from
%them that only contains values in the cloudy layers - am3_lwc_cloud
make_indices_for_continuous_k

%am3_lwc_cloud is in g/kg at this stage
am3_mean_lwc_cloud = meanNoNan(am3_lwc_cloud,1); %mean over each vertical dim now gives the mean lwc
[me,Nme]=meanNoNan(1e3*am3_lwc_cloud.*permute(am3_rho,[2 1 3 4]).*-dz,1); %mean over each vertical dim now gives the mean lwc
am3_LWC_max = 1e3*max(am3_lwc_cloud.*permute(am3_rho,[2 1 3 4]),[],1); %max LWC within the designated cloudy layer
am3_LWC_max = squeeze(am3_LWC_max);
am3_lwp_cloud = me.*Nme; %equivalent to sum(xxx)
%LWP is in g/m2

[am3_Nd_meanliq,Nlevs_meanNd] = meanNoNan(am3_Nd_cloud,1);
am3_Nd_max = max(am3_Nd_cloud,[],1); %max Nd within the designated cloudy layer
ilevs = find(Nlevs_meanNd<thresh_Nlevs);
am3_Nd_max = squeeze(am3_Nd_max);
am3_Nd_meanliq(ilevs)=NaN;
am3_Nd_max(ilevs)=NaN;
am3_LWC_max(ilevs)=NaN;

am3_CF_max = squeeze(max(am3_CF_cloud,[],1));


%IMPORTANT! - remove values when the max or min was NaN (i.e. all the values are NaN)
%- since an index will be reported, but it is meaningless
am3_liq_maxliq(isnan(max_liq))=NaN;
%am3_Nd_maxliq(isnan(max_liq))=NaN;
am3_Nd_maxliq(isnan(am3_CTP))=NaN;

am3_CF_maxliq(isnan(max_liq))=NaN;
am3_CTH(isnan(am3_CTP))=NaN;
am3_CBH(isnan(am3_CTP))=NaN;
am3_CBP(iremove2)=NaN;



%permute them back to the normal way around
am3_liq2 = permute(am3_liq2,[2 1 3 4]);
am3_drop2 = permute(am3_drop2,[2 1 3 4]);
am3_cf2 = permute(am3_cf,[2 1 3 4]);
h_half2 = permute(h_half2,[2 1 3 4]);
am3_pfull2 = permute(am3_pfull2,[2 1 3 4]);

disp('Done am3_process')

%new function 
%        [ilat,ilon,dmins]=am3_map_lat_lon_points(am3_lat,am3_lon,lat_am3,lon_am3);

%would still need to calculate Plat2 and Plon2 below and feed these into
%the above

% lat_min = min(am3_lat,[],1);
% lat_max = max(am3_lat,[],1);
% 
% Plat2=[ceil(min(lat_max)):1:floor(max(lat_min))];
% Plon2=[ceil(minALL(am3_lon)):1:floor(maxALL(am3_lon))];
% 
% clear am32_Nd am32_CF dmins
% max_ov=-1;
% for i=1:length(Plat2)
%     for j=1:length(Plon2)
%         dist = distlatlon(Plat2(i),Plon2(j),am3_lat,am3_lon);
%         [minval,idist]=minALL(dist);
%         lats(i,j) = am3_lat(idist(1),idist(2));
%         lons(i,j) = am3_lon(idist(1),idist(2));
%         
%         am32_Nd(i,j) = am3_Nd_maxliq(idist(1),idist(2));
%         am32_CF(i,j) = max(am3_cf(:,idist(1),idist(2)),[],1);        
% 
%         
%         dmin = minALL(dist);
%         if dmin>max_ov
%             max_ov=dmin;
%             imax=i;
%             jmax=j;
%         end
%         
%         dmins(i,j)=dmin;
% 
% 
%     end
% end
% 
% [Plon,Plat]=meshgrid(Plon2,Plat2);













