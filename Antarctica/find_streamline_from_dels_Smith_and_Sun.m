function [z]=find_streamline_from_dels_Smith_and_Sun(z0,del_all,z_d_arr)
%function [z]=find_streamline_from_dels_Smith_and_Sun(z0,del_all,z_d_arr)
%finds the height of a streamline starting at z_hat=z0
%based on the del_hat field del_all and the height_hat field z_d_all

for i=1:size(del_all,1)
    z_a=z0+del_all(i,:); %height that the lower A streamline should be at
    [temp,iz]=min(abs(z_a-z_d_arr(i,:)));
    z(i)=z_d_arr(i,iz);
end