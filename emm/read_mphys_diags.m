function [mphys]=read_mphys_diags(emmdir2)

mp=dlmread([emmdir2 'mphys_diags'],' ');

a=find(mp(:,1)==0.5);
s1=length(a);
b=size(mp,1);
s2=floor(b/s1);

mphys=mp(1:s1*s2,1);
mphys=reshape(mphys,[s1 s2]); %updraught
