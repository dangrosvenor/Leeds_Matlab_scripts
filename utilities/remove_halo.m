function Y=remove_halo(X,halo_size)



sX2 = size(X,2);

Y = X(:,halo_size+1:sX2-halo_size);

