h=fspecial('average',[30 30]);
nonan=vals;
b=isnan(vals);
bb=find(b==1);
nonan(bb)=0; %replace nans with zero

smoo=filter2(h,nonan); %smooth out data with no nans

b=find(smoo<10); %data to be ignored

v1=size(vals,1);
v2=size(vals,2);
x=repmat([1:v1]',[1 v2]);
y=repmat([1:v2],[v1 1]);

[th,r]=cart2pol(x-241,y-241);
br=find(r>50);

newim=vals;
newim(b)=0;
newim(br)=nonan(br);


%figure;pcolor(newim);caxis([0 60]);colorbar;shading flat;axis ij