function V2=wrap3d(V2)
%program to wrap 3d fields to make the cloud in the centre

[ix,iy,iz]=size(V2);
ixh=ceil(ix/2);
iyh=ceil(iy/2);

%iyh=60;

V(1:ix-ixh,:,:)=V2(ixh+1:ix,:,:);
V(ix-ixh+1:ix,:) = V2(1:ixh,:);
%V=V';

V2(:,1:iy-iyh,:,:)=V(:,iyh+1:iy,:);
V2(:,iy-iyh+1:iy,:) = V(:,1:iyh,:);
%V=V';
