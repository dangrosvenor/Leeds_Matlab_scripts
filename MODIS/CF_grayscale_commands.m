
%set the screening for ihtot to NOT include CF
%we make up the full screening set later
ihtot_other = ihtot; %store the ihtot values for the things
%other than CF

if icolormap_cf_grey==1
    %points where just the CF is bad
    ihtot_cf = find( ~ ( cf_time3>=thresh_CF) );

    %the points where either CF or the others are bad (full screened set)
    ihtot = union(ihtot_cf,ihtot); %union combines with no repetitions
end

%remove points we defo don't want for Nd calc (full
%screening)
dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

%do the mean over all the required times and for
%both satellites if selected (Nd mean)
dat_modis_temp = dat_modis(:,:,time_inds_average);
if length(size(dat_modis_temp))==3
[P,Npoints] = meanNoNan(dat_modis_temp,3); %time mean (all times)
else
    P = dat_modis_temp;
    Npoints = ones*size(dat_modis_temp);
end
%and now do CF screening stuff - In case where one satellite
%gives a high CF and so a valid Nd, but the
%other does not we want to use the Nd from the good satellite
%we only do CF if both satellites are bad

if icolormap_cf_grey==1
    %make a NaN array - this will only include the low
    %CF points
    dat_modis2 = NaN*ones(size(dat_modis));
    %scale the cf to match the CF gray colorbar added
    %on the end
    dat_modis2(ihtot_cf) = 1e12+1e21*cf_time3(ihtot_cf);
    %remove the points that are bad for other reasons
    %e.g. don't want to include the CF for a bad sensor
    %ZA
    dat_modis2(ihtot_other)=NaN;

    %new P average - only for low CF points, that are not screened
    %for other reasons. everything else is NaN & will be ignored
    [P2,Npoints2] = meanNoNan(dat_modis2(:,:,time_inds_average),3); %time mean (all times)

    %we only want to replace points that were screened
    %before - i.e. not good points where we have
    %droplet numbers from one or both satellites
    inanP = find(isnan(P)==1);
    P(inanP)=P2(inanP);

end


