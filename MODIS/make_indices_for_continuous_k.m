dK = ( iCBP(:)-iCTP(:) + 1);  %since the indices are backwards low means higher altitude
i0 = find(iCBP>1); %if iCBP==1 then this will be a non-cloudy column
NK = sum(dK(i0));

Ks = NaN*ones([1 NK]);
Is = NaN*ones([1 NK]);

loc=1;
for i=1:length(i0)
    K = iCTP(i0(i)):iCBP(i0(i));
    loc2=loc+length(K)-1;
    Ks(loc:loc2) = K;
    Is(loc:loc2) = i0(i);
    loc=loc2+1;
end

%convert our K and I index into linear indices in e.g. gcm_lwc_minthreshCF2
icloud_lin = sub2ind([s_gcm_low(1) s_gcm_low(2)*s_gcm_low(3)*s_gcm_low(4)],Ks(:),Is(:));

gcm_lwc_cloud = NaN*ones(size(gcm_lwc_minthreshCF2));
gcm_lwc_cloud(icloud_lin) = gcm_lwc_minthreshCF2(icloud_lin);

gcm_Nd_cloud = NaN*ones(size(gcm_lwc_minthreshCF2));
%gcm_Nd_cloud(icloud_lin) = gcm_drop2(icloud_lin);
gcm_Nd_cloud(icloud_lin) = gcm_drop_ice_screen(icloud_lin);

gcm_REFFL_cloud = NaN*ones(size(gcm_lwc_minthreshCF2));
gcm_REFFL_cloud(icloud_lin) = gcm_REFFL_ice_screen(icloud_lin);


gcm_CF_cloud = NaN*ones(size(gcm_lwc_minthreshCF2));
gcm_CF_cloud(icloud_lin) = gcm_CF_screened(icloud_lin);

if exist('h_half2')
    dz = diff(h_half2,1,1); %thickness of the model layers
end

