function Y=add_halo(X,halo_size)

s2=size(X,2);
Y=NaN*ones(size(X,1),size(X,2)+halo_size*2);


sY2=size(Y,2);
Y(:,1:halo_size) = X(:,s2-halo_size+1:s2);
Y(:,sY2-halo_size+1:sY2) = X(:,1:halo_size);
Y(:,halo_size+1:halo_size+s2) = X;