%plots 3d iso images liquid=blue ice=red graupel=yellow snow=grey

cen=0; %cen=1 for centring
xoff=0; %was=9000
yoff=0;

clear x y z v

figure;

[a b c]=size(ThreeD.Q(:,:,:,6));

tol=50;
tt=[1 1];
stt=size(tt);
while stt(2)>1 
    tt=find(abs(Grid.X1-Grid.X1(1)-xoff)<tol);
    stt=size(tt);
    tol=tol+50;
end
ttx=tt;

tt=[1 1];
stt=size(tt);
while stt(2)>1 
    tt=find(abs(Grid.Y1-Grid.Y1(1)-yoff)<tol);
    stt=size(tt);
    tol=tol+50;
end
tty=tt;



max(ThreeD.W(:,:,:),[],3);
[a3 b3]=max(max(ans,[],1));

if b3(1)<a/2
    xbmax=b3(1)+ttx;
else
    xbmax=b3(1)-ttx;
end

max(ThreeD.Q(:,:,:,6),[],3);
[a3 b3]=max(max(ans,[],2));

if b3(1)<b/2
    ybmax=b3(1)+tty;
else
    ybmax=b3(1)-tty;
end





xs=xbmax-(0.5*(a));
ys=ybmax-(0.5*(b));




z=Grid.Z;

if xs<0
    xr=-xs;
    xs=(0.5*(a))+xbmax;
else
        xr=(1.5*(a))-xbmax;
end
    
if ys<0
    'yesiso'
    yr=-ys;
    ys=(0.5*(b))+ybmax;
else
    yr=(1.5*(b))-ybmax;
end   

%x(1:xr-1)=Grid.X1(xs+1:a-1);
%x(xr:a-2)=Grid.X1(2:xs);


% x=Grid.X1(2:b-1);
% y(1:yr-1)=Grid.Y1(ys+1:b-1);
% y(yr:b-2)=Grid.Y1(2:ys);
% xs=1;
% xr=a-1;
% ys=1;
% yr=b-1;
% 
% x(1:xs-1)=Grid.X1(xr+1:a-1);
% x(xs:a-2)=Grid.X1(2:xr);
% 
% 
% %x=Grid.X1(2:b-1);
% y(1:ys-1)=Grid.Y1(yr+1:b-1);
% y(ys:b-2)=Grid.Y1(2:yr);



% ThreeD2.Q(1:xr,:,:,:)=ThreeD.Q(xs+1:a,:,:,:);
% ThreeD2.Q(xr+1:a,:,:,:)=ThreeD.Q(1:xs,:,:,:);
% 
% ThreeD2.Q(:,1:yr,:,:)=ThreeD.Q(:,ys+1:b,:,:);
% ThreeD2.Q(:,yr+1:b,:,:)=ThreeD.Q(:,1:ys,:,:);

x=Grid.X1;
y=Grid.Y1;


if cen==1
v2(1:xr,:,:)=ThreeD.Q(xs+1:a,:,:,2); %Liquid = blue
v2(xr+1:a,:,:)=ThreeD.Q(1:xs,:,:,2);

v(:,1:yr,:)=v2(:,ys+1:b,:);
v(:,yr+1:b,:)=v2(:,1:ys,:);
else
   v=ThreeD.Q(1:a,1:b,:,2);
end

colour='blue';
plotiso3d;
view(3);
%axis tight;
camlight;

if cen==1
v2(1:xr,:,:)=ThreeD.Q(xs+1:a,:,:,5); %graupel = yellow
v2(xr+1:a,:,:)=ThreeD.Q(1:xs,:,:,5);

v(:,1:yr,:)=v2(:,ys+1:b,:);
v(:,yr+1:b,:)=v2(:,1:ys,:);
else
v=ThreeD.Q(1:a,1:b,:,5);
end

colour='yellow';
plotiso3d;

if cen==1
v2(1:xr,:,:)=ThreeD.Q(xs+1:a,:,:,4); %snow = grey
v2(xr+1:a,:,:)=ThreeD.Q(1:xs,:,:,4);

v(:,1:yr,:)=v2(:,ys+1:b,:);
v(:,yr+1:b,:)=v2(:,1:ys,:);
else
v=ThreeD.Q(1:a,1:b,:,4);
end

colour=[0.9 0.9 0.9];
plotiso3d;

if cen==1
v2(1:xr,:,:)=ThreeD.Q(xs+1:a,:,:,6); %ice = red
v2(xr+1:a,:,:)=ThreeD.Q(1:xs,:,:,6);

v(:,1:yr,:)=v2(:,ys+1:b,:);
v(:,yr+1:b,:)=v2(:,1:ys,:);
else
v=ThreeD.Q(1:a,1:b,:,6);
end

colour='red';
plotiso3d;

axis([Grid.X1(1) Grid.X1(end) Grid.Y1(1) Grid.Y1(end) Grid.Z(1) Grid.Z(end)]);


% v=ThreeD.Q(:,:,:,2);
% p=patch(isosurface(x,y,z,v));
% 
% lighting gouraud;
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
% 
% hold on;
% 
% v=ThreeD.Q(:,:,:,5);
% p=patch(isosurface(x,y,z,v));
% set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
% lighting gouraud;
% 
% v=ThreeD.Q(:,:,:,6);
% p=patch(isosurface(x,y,z,v));
% set(p, 'FaceColor', 'yellow', 'EdgeColor', 'none');
% lighting gouraud;

