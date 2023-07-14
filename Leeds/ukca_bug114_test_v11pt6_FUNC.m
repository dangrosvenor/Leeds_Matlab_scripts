function [zNd,Nd_prof,i,j] = ukca_bug114_test_v11pt6_FUNC(um_case,ioverride_ij,i,j)

if ~exist('ioverride_ij')
   ioverride_ij=0; 
end

var_UM = 'Nd_3d_38479';
%um_case='u-bt398'; 
dirUM = ['/home/disk/eos15/d.grosvenor/UM/12Nov2008_Boutle/' um_case];

%filename = [dirUM '/' var_UM '/umnsaa_pb012_201605101300_' var_UM '_saved.nc'];
filename = [dirUM '/' var_UM '/merged.nc'];

nc = netcdf(filename);
Nd = nc{var_UM};
zNd = nc{'height'}(:);

z=3000;

[tempval,iz]=min(abs(z-zNd));
it=size(Nd,1);

if ioverride_ij==0
    
    Nd_2d = Nd(it,iz,:,:);
    [tempval,ind]=max(Nd_2d(:),[],1);
    
    [i,j] = ind2sub(size(Nd_2d),ind);
end
Nd_prof = Nd(it,:,i,j);


nc = close(nc);